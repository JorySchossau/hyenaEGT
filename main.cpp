//
//  main.cpp
//  PTG
//
//  Created by Arend Hintze on 4/11/14.
//  Copyright (c) 2014 Arend Hintze. All rights reserved.
//  Modified by Jory Schossau 2015-2016
//
#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <map>
#include <cmath> // for pow
#include <ctime> // for time
#include <iostream> // for console display
#include <iomanip> // for formatting decimals
#include <curses.h>
#include <functional> // for ref
#include <random>
#include "params/params.h"

#ifdef _WIN32
	#include <process.h> // process ID
   #include <direct.h> // filesystem access
#else
	#include <unistd.h> // process ID
   #include <sys/stat.h> // filesystem access
#endif

using namespace Params;

bool fileExists(const string& name) {
   struct stat buffer;
   return( stat( name.c_str(), &buffer) == 0); // 0 if no problem and file is existing
}

//void computeLOD(FILE *f,FILE *g, tAgent *agent,tGame *game);
const char* cstr(string s) { return s.c_str(); }
double randDouble() { return ((double)rand() / (double)RAND_MAX); }

#define C 1
#define D 0
#define M 3
#define I 2

#define genes 2
// xDim, yDim assumed powers of 2
#define xDim 32
#define yDim 32
#define possibleMoves 4

// neighbors is assumed power of 2
#define neighbors 4
int xm[neighbors+1]={0,1,0,-1,0};
int ym[neighbors+1]={-1,0,1,0,0};
//int xm[neighbors+1]={0,1,0,-1, 1,1,-1,-1, 0,2,0,-2, 2,2,-2,-2, 0};
//int ym[neighbors+1]={-1,0,1,0, 1,-1,-1,1, -2,0,2,0, 2,-2,-2,2, 0};
vector< float > ppc; // Pc & Pd inheritances: 0.0 <= p <= 1.0 and is constant value otherwise -1 implies usual mutation
vector< float > fcdmi; // fraction of initial C,D,M,I players

namespace g {
	string experimentID;
	int replicateID;

	double mutationRate;
	double replacementRate=0.02;

	double beta; // cost of being punished
	double gamma; // cost of punishing
	float rMultiplier; // synergy
	double zeta; // threshold limit
	int radius=32; // radius of play. for world xDim of 32, radius 32 is well-mixed
	int updates=42; // usually in thousands
	bool savePlays; // Masoud's request for information calculation
	bool deterministic; // makes game deterministic
	bool thresholdPayoff;
	int update=0;
	int clanSize=5; // defunct
	float cost; // cost of cooperation
}

class tPlayer{
public:
   double probs[genes];
   unsigned char action;
   tPlayer *ancestor;
   int nrPointingAtMe;
   float score;
   int born;
   vector<char> plays;
   
   tPlayer();
   ~tPlayer();
   void inherit(tPlayer *from);
   unsigned char move(void);
   double inheritGene(int w);
   void exportLOD(FILE *f);
   char* getPlayHistoryAsString();
   const char playToChar[4] = {'c','d','m','i'}; // lookup table for unsigned char (int) to readable char
};

struct PopulationFrequencies {
   float c,d,m,i;
   void normalize() {
      float total = c+d+m+i;
      c/=total;
      d/=total;
      m/=total;
      i/=total;
      return;
   }
	void round() {
		if (c < 1) c = 0;
		if (d < 1) d = 0;
		if (m < 1) m = 0;
		if (i < 1) i = 0;
		return;
	}
};

namespace Threshold {
	enum type {Stepwise, StepwiseHalf, Sigmoid}; // other types StepwiseHalf, Sigmoid
}

inline float threshold(const float val, Threshold::type type = Threshold::Stepwise) {
	float result(0.0);
	switch(type) {
		case Threshold::Stepwise:
			result = (val < 0.0) ? 0.0 : 1.0;
			break;
		case Threshold::StepwiseHalf:
			result = (val < 0.0) ? 0.0 : 1.0;
			if (val == 0.0) result = 0.5;
			break;
		case Threshold::Sigmoid:
			result = 1.0 / ( 1.0 + exp( -val * 100.0 ) );
			break;
	}
	return result;
}


int main (int argc, char* argv[]) {
   int x,y,z,i,tx,ty;
   double maxFit,minFit;
	int maxFitX(0), maxFitY(0);
   vector< vector<tPlayer*> > player;
   vector< PopulationFrequencies* > census; /// record of frequencies over time
	PopulationFrequencies* frequencies;
   int lxm[neighbors],lym[neighbors];
   srand(getpid());
	FILE *fileSRAND;
	//fileSRAND=fopen(cstr("srand.txt"),"w+a");
	//fprintf(fileSRAND, "%d\n", getpid());
	//fclose(fileSRAND);
   bool debug;
	double pool(0.0);
	unsigned char N[neighbors], done[neighbors+1];

   bool showhelp;
	bool cursesMode; // for visual display of game
   /// Filenames
#pragma optimize off
   string filenameLOD;
#pragma optimize on
   string filenameGenome;
   string filenamePOP;
   string filenameEnd;
   FILE *fileLOD;
   FILE *fileGenome;
   FILE *filePOP;
   FILE *fileEnd;
	float param(0.0);
   /// other parameters

   
   /// required
	addp(TYPE::STRING, &g::experimentID, "--experiment", "unique identifier for this experiment, shared by all replicates.");
	addp(TYPE::INT, &g::replicateID, "--replicate", "unique number to identify this replicate in this experiment.");
   /// not required
	addp(TYPE::BOOL, &showhelp);
   addp(TYPE::BOOL, &debug, "false", false, "--debug", "enables verbose printing.");
   addp(TYPE::BOOL, &g::savePlays, "false", false, "--savePlays", "Exports play history in line of descent output.");
   addp(TYPE::BOOL, &g::deterministic, "false", false, "--deterministic", "Forces a deterministic strategy simulation (random initial pop).");
   addp(TYPE::DOUBLE, &g::beta, "0.0", false, "--beta", "Beta parameter: effect of punishment.");
   addp(TYPE::DOUBLE, &g::gamma, "0.0", false, "--gamma", "Gamma parameter: cost of punishment.");
   addp(TYPE::FLOAT, &g::rMultiplier, "5", false, "--r", "R multiplier parameter.");
   addp(TYPE::DOUBLE, &g::zeta, "0.0", false, "--zeta", "Payoff zeta for group success.");
	addp(TYPE::FLOAT, &g::cost, "1.0", false, "--cost", "Cost of cooperation.");
   addp(TYPE::INT, &g::radius, "32", false, "--radius", "Length of one side of square group region, 32 is K=1024 which implies well-mixed since popoulation is 1024.");
	addp(TYPE::STRING, &filenameLOD, "none", false, "--lod", "filename to save Line of Descent.");
	addp(TYPE::STRING, &filenameGenome, "none", false, "--genome", "filename to save LCA genome.");
	addp(TYPE::STRING, &filenamePOP, "none", false, "--pop", "filename to save population frequencies.");
   addp(TYPE::STRING, &filenameEnd, "none", false, "--end", "filename to save best genome and pop frequencies.");
	addp(TYPE::INT, &g::updates, "100000", false, "--updates", "number of updates to simulate (not generations).");
	addp(TYPE::INT, &g::clanSize, "5", false, "--clanSize", "size of one clan (group) within which are C & D strategies.");
   addp(TYPE::FLOAT, &ppc, 2, "-1.0", false, "--ppc", "The probabilities for punishment and cooperation (0-1)");
	addp(TYPE::FLOAT, &fcdmi, 4, "-1.0", false, "--fcdmi", "The initial population frequencies for C D M I (0-1)");
   addp(TYPE::DOUBLE, &g::mutationRate, "0.02", false, "--mu", "mutation rate (per site)");
	addp(TYPE::BOOL, &g::thresholdPayoff, "false", false, "--thresholdPayoff", "Uses the thresholding fn based on Nc to determine if the group receives any payoff.");
	addp(TYPE::BOOL, &cursesMode, "false", false, "--curses", "Uses ncurses to display the game.");

	argparse(argv);
	if (showhelp) {
		std::cout << argdetails() << std::endl;
		std::cout << "Example minimal invocation:" << std::endl;
		std::cout << argv[0] << " --experiment=deterministic --replicate=1 --lod=lineOfDescent.lod --genome=genome.gen" << std::endl;
		std::cout << "or" << std::endl;
		std::cout << argv[0] << " --experiment deterministic --replicate 1 --lod lineOfDescent.lod --genome genome.gen" << std::endl;
		std::cout << "or" << std::endl;
		std::cout << argv[0] << " --experiment trial --replicate 1 --end out.end --deterministic --beta 0 --gamma 0 --fcdmi 0.25 0.25 0.25 0.25 --debug --updates 10000 --zeta 0 --r 3" << std::endl;
		std::cout << std::endl;
		std::cout << string("There are ") + to_string(neighbors) + string(" neighbors defined in this build.") << std::endl;
		std::cout << std::endl;
		exit(0);
	}

	if (cursesMode) {
		// initialize the screen
		initscr();
		curs_set(0);
		if (not has_colors() ) std::cout << "WARNING: colors unsupported" << std::endl;
		else start_color();
		init_pair(1, COLOR_GREEN, COLOR_GREEN);
		init_pair(2, COLOR_BLUE, COLOR_BLUE);
		init_pair(3, COLOR_WHITE, COLOR_GREEN);
		init_pair(4, COLOR_WHITE, COLOR_BLUE);
		init_pair(5, COLOR_BLACK, COLOR_GREEN);
		init_pair(6, COLOR_BLACK, COLOR_BLUE);
	}

	/// initialize player matrix
   player.resize(xDim);
   for(x=0;x<xDim;x++) {
      player[x].resize(yDim);
      for(y=0;y<yDim;y++) {
         player[x][y]=new tPlayer();
		}
   }
   for(g::update=1;g::update<g::updates;g::update++){
      frequencies = new PopulationFrequencies;
      frequencies->c=0; frequencies->d=0; frequencies->m=0; frequencies->i=0; 
		/// single-threaded version
		param = 0.0;
      for(x=0;x<xDim;x++) {
         for(y=0;y<yDim;y++) {
				////if ((g::update<2048) && ((x|y)==0)) player[x][y]->probs[1]=1.0;
				//static int n(0); // number of hunters clan produces
				//n = player[x][y]->probs[1]*g::clanSize; // Use that genome's probability to determine number of cooperators (hunters)
				//if ( randDouble() < threshold(n-g::zeta, Threshold::Stepwise) ) { // determine group's score
				//	player[x][y]->score += n*((static_cast<double>(g::rMultiplier)/static_cast<double>(g::clanSize)) - 1.0) + (g::clanSize-n)*(static_cast<double>(g::rMultiplier)/static_cast<double>(g::clanSize)); // C + D payoffs when above threshold
				//	frequencies->c+=1;
				//} else {
				//	player[x][y]->score += n*(-1) + (g::clanSize-n)*(0); // C + D payoffs when below threshold
				//	frequencies->d+=1;
				//}
				////if (player[x][y]->score < 0.0) player[x][y]->score = 0.0;
				//param+=player[x][y]->score;
				////if (player[x][y]->score > param) param=player[x][y]->score;
				//
				//if ((g::update&127)==127) {
				//	if (cursesMode) {
				//		if (threshold(n-g::zeta, Threshold::Stepwise)) {
				//			mvaddch(y,x*2,'.'|COLOR_PAIR(1));
				//			addch('.'|COLOR_PAIR(1));
				//		} else {
				//			mvaddch(y,x*2,'.'|COLOR_PAIR(2));
				//			addch('.'|COLOR_PAIR(2));
				//		}
				//	}
				//}
				for(i=0;i<neighbors;i++){ /// create a group to play games
					xm[i] = (rand()&(g::radius-1)) * ((rand()&2)-1);
					ym[i] = (rand()&(g::radius-1)) * ((rand()&2)-1);
				}
            for(z=0;z<possibleMoves;z++) N[z]=0; /// reset group's moves
            for(z=0;z<neighbors+1;z++){ /// everyone makes a play, even 'this' player (last neighbor)
               done[z]=player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->move();
               player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->plays.push_back(done[z]); // record play
               N[done[z]]++;
            }
            frequencies->c+=N[C]; frequencies->d+=N[D]; frequencies->m+=N[M]; frequencies->i+=N[I];
				if (g::thresholdPayoff) {
					if ( randDouble() < threshold(N[C]+N[M]-g::zeta, Threshold::Stepwise) ) pool = 1.0;
					else pool = 0.0;
				} else {
					pool = N[C]+N[M];
				}
            for(z=0;z<neighbors+1;z++) {
               switch(done[z]){
                  case C: //C ooperator
                     player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score+=( (g::rMultiplier)*((pool)/((double)neighbors+1.0))-g::cost);
                     break;
                  case D: //D efector
                     player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score+=( (g::rMultiplier)*((pool)/((double)neighbors+1.0))-(g::beta*((double)N[M]+(double)N[I])/((double)neighbors)));
                     break;
                  case M://M oralist
                     player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score+=( (g::rMultiplier)*((pool)/((double)neighbors+1.0))-g::cost-(g::gamma*((double)N[D]+(double)N[I])/((double)neighbors)));
                     break;
                  case I://I moralist
                     player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score+=( (g::rMultiplier)*((pool)/((double)neighbors+1.0))-(g::beta*((double)N[M]+(double)N[I]-(double)1.0)/((double)neighbors))-(g::gamma*((double)N[D]+(double)N[I])/((double)neighbors)));
                     break;
               }
            }
			}
		}
		if (cursesMode && ((g::update&4095)==4095)) {
			refresh();
		}
		for(z=0;z<int(xDim*yDim*g::replacementRate);z++){
			do{ /// pick player to die
				x=rand()&(xDim-1);
				y=rand()&(yDim-1);
			} while(player[x][y]->born==g::update);

			maxFit=player[0][0]->score;
			minFit=maxFit;
			for(x=0;x<xDim;x++) {
				for(y=0;y<yDim;y++) {
					maxFit=fmax(maxFit,player[0][0]->score);
					minFit=fmin(minFit,player[0][0]->score);
				}
			}
			//for(tx=0;tx<g::radius;tx++) { /// determine fitness range for pool of reproduction candidates
			//	for(ty=0;ty<g::radius;ty++) {
			//		if(player[(x+tx)&(xDim-1)][(y+ty)&(yDim-1)]->score>maxFit) /// ++
			//			maxFit=player[(x+tx)&(xDim-1)][(y+ty)&(yDim-1)]->score;
			//		if(player[(x-tx)&(xDim-1)][(y-ty)&(yDim-1)]->score<minFit) /// --
			//			maxFit=player[(x-tx)&(xDim-1)][(y-ty)&(yDim-1)]->score;
			//		if(player[(x+tx)&(xDim-1)][(y-ty)&(yDim-1)]->score>maxFit) /// +-
			//			maxFit=player[(x+tx)&(xDim-1)][(y-ty)&(yDim-1)]->score;
			//		if(player[(x-tx)&(xDim-1)][(y+ty)&(yDim-1)]->score<minFit) /// -+
			//			maxFit=player[(x-tx)&(xDim-1)][(y+ty)&(yDim-1)]->score;

			//		if(player[(x+tx)&(xDim-1)][(y+ty)&(yDim-1)]->score>maxFit) /// ++
			//			minFit=player[(x+tx)&(xDim-1)][(y+ty)&(yDim-1)]->score;
			//		if(player[(x-tx)&(xDim-1)][(y-ty)&(yDim-1)]->score<minFit) /// --
			//			minFit=player[(x-tx)&(xDim-1)][(y-ty)&(yDim-1)]->score;
			//		if(player[(x+tx)&(xDim-1)][(y-ty)&(yDim-1)]->score>maxFit) /// +-
			//			minFit=player[(x+tx)&(xDim-1)][(y-ty)&(yDim-1)]->score;
			//		if(player[(x-tx)&(xDim-1)][(y+ty)&(yDim-1)]->score<minFit) /// -+
			//			minFit=player[(x-tx)&(xDim-1)][(y+ty)&(yDim-1)]->score;
			//	}
			//}
			if ((maxFit-minFit) <= 0.0) { /// Fitness-proportionally select a player to reproduce
				do{ /// choose reproducer when no good choices
					lxm[0] = (rand()&(g::radius-1)) * ((rand()&2)-1);
					lym[0] = (rand()&(g::radius-1)) * ((rand()&2)-1);
				} while(((lxm[0]==0)&&(lym[0]==0)) || (player[(x+lxm[0])&(xDim-1)][(y+lym[0])&(yDim-1)]->born==g::update));
			} else { /// maxFit > 0.0
				do{ /// choose reproducer when there are fitness differences
					lxm[0] = (rand()&(g::radius-1)) * ((rand()&2)-1);
					lym[0] = (rand()&(g::radius-1)) * ((rand()&2)-1);
				} while( (randDouble()>((player[(x+lxm[0])&(xDim-1)][(y+lym[0])&(yDim-1)]->score-minFit)/(maxFit-minFit))) || ((lxm[0]==0)&&(lym[0]==0)) || (player[(x+lxm[0])&(xDim-1)][(y+lym[0])&(yDim-1)]->born==g::update) );
			}
         
			player[x][y]->nrPointingAtMe--;
			if(player[x][y]->nrPointingAtMe==0)
				delete player[x][y];
			player[x][y]=new tPlayer;
			player[x][y]->inherit(player[(x+lxm[0])&(xDim-1)][(y+lym[0])&(yDim-1)]);
      }
      frequencies->normalize();
      census.push_back(frequencies);
		if (!cursesMode) {
			if (debug && (g::update&511) == 511) {
				tPlayer* ptr = player[0][0];
				for (i=0; i<500; i++) {
					if (ptr->ancestor == nullptr) break;
					ptr = ptr->ancestor;
				}
				std::cout << std::fixed;
				std::cout<<std::setprecision(2)<<frequencies->c<<" "<<frequencies->d<<" "<<frequencies->m<<" "<<frequencies->i<<"\t"<<ptr->probs[0]<<"\t"<<ptr->probs[1]<<"\t"<<param<<"\t"<<param/(xDim*yDim)<<std::endl;
			}
			//if (frequencies->c == 1.0 || frequencies->d == 1.0 || frequencies->m == 1.0 || frequencies->i == 1.0) {
			//	if (debug) printf( "%f %f %f %f\n",frequencies->c, frequencies->d, frequencies->m, frequencies->i );
			//	break;
			//}
		}
   }
   if (filenameLOD != "none") {
      fileLOD=fopen(cstr(filenameLOD),"w+t");
      if (g::savePlays)
         fprintf( fileLOD, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
               "experiment", "replicate", "born",
               "r","zeta","beta","gamma",
               "p0","p1",
               "pc","pd","pm","pi",
               "plays" );
      else
         fprintf( fileLOD, "%s %s %s %s %s %s %s %s %s %s %s %s %s\n",
               "experiment", "replicate", "born",
               "r","zeta","beta","gamma",
               "p0","p1",
               "pc","pd","pm","pi" );
      player[0][0]->exportLOD(fileLOD);
      fclose(fileLOD);
   }
   if (filenamePOP != "none") {
      filePOP=fopen(cstr(filenamePOP), "w+t");
      fprintf( filePOP, "%s %s %s %s %s\n","update","c","d","m","i" );
      for (int pop_i=0; pop_i<census.size(); ++pop_i) {
         fprintf( filePOP, "%d %f %f %f %f\n",pop_i,census[pop_i]->c,census[pop_i]->d,census[pop_i]->m,census[pop_i]->i );
      }
      fclose(filePOP);
   }
   if (filenameEnd != "none") {
		if (player[0][0]->ancestor!=nullptr && player[0][0]->ancestor->ancestor!=nullptr && player[0][0]->ancestor->ancestor->ancestor!=nullptr) {
			/// write endpoint data to file, appending data,
			/// otherwise creating (with csv header) if necessary
			/// before appending
			bool fileEndExists = false;
			if ( fileExists(filenameEnd) ) fileEndExists = true;
			fileEnd=fopen(cstr(filenameEnd), "a+t");
			if (fileEndExists == false) {
				fprintf( fileEnd, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
						"experiment", "replicate",
						"r","zeta","beta","gamma",
						"c","d","m","i",
						"p0","p1",
						"pc","pd","pm","pi",
						"fc","fd","fm","fi" );
			}
			tPlayer* LCA = player[0][0]->ancestor->ancestor->ancestor;
			for (i=0; i<500; i++) { // get an older LCA
				if (LCA->ancestor == nullptr) break;
				LCA = LCA->ancestor;
			}
			fprintf(fileEnd, "%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
					cstr(g::experimentID), g::replicateID,
					g::rMultiplier, g::zeta, g::beta, g::gamma,
					frequencies->c, frequencies->d, frequencies->m, frequencies->i,
					LCA->probs[0], LCA->probs[1],
					(1.0-LCA->probs[0])*LCA->probs[1],
					//cgene_value,
					(1.0-LCA->probs[0])*(1.0-LCA->probs[1]),
					LCA->probs[0]*LCA->probs[1],
					LCA->probs[0]*(1.0-LCA->probs[1]),
					fcdmi[0], fcdmi[1], fcdmi[2], fcdmi[3]);
			fclose(fileEnd);
		}
   }
	if (cursesMode) endwin();
   return 0;
}

tPlayer::tPlayer() {
   int i;
	if (g::deterministic) {
		double a(0.25),b(0.25),c(0.25),d(0.25);
		if (fcdmi[0] > 0.0) {
			a = fcdmi[0];
			b = fcdmi[1];
			c = fcdmi[2];
			d = fcdmi[3];
		} else {
			a = randDouble();
			b = randDouble();
			c = randDouble();
			d = randDouble();
			double total=a+b+c+d;
			a/=total;
			b/=total;
			c/=total;
			d/=total;
		}
		double rnd = randDouble();
		if (rnd <= a) { probs[0]=0.0; probs[1]=1.0; }
		else if (rnd <= (a+b)) { probs[0]=0.0; probs[1]=0.0; }
		else if (rnd <= (a+b+C)) { probs[0]=1.0; probs[1]=1.0; }
		else { probs[0]=1.0; probs[1]=0.0; }
	} else {
		for(i=0;i<genes;i++) {
			if (ppc[i] < 0.0) {
				probs[i]=randDouble();
			} else {
				probs[i]=ppc[i];
			}
		}
	}
	if (g::deterministic) { // if exact player types specified
		action=0;
		for(i=0;i<genes;i++)
			if(randDouble()<probs[i])
				action=(action<<1)+1;
			else
				action=action<<1;
   }
   score=0;
   ancestor=NULL;
   nrPointingAtMe=1;
   born=g::update;
}

tPlayer::~tPlayer() {
   if(ancestor!=NULL) {
      ancestor->nrPointingAtMe--;
      if(ancestor->nrPointingAtMe==0)
         delete ancestor;
   }
}

unsigned char tPlayer::move(void) {
   unsigned char r=0;
   int i;
   if (g::deterministic) {
      r=action;
   } else {
   for(i=0;i<genes;i++)
      if(randDouble()<probs[i])
         r=(r<<1)+1;
      else
         r=r<<1;
	}
   return r;
}
//0.0, 0.0 => 11 => M
//0.0, 1.0 => 10 => I
//1.0, 0.0 => 01 => C
//1.0, 1.0 => 00 => D
void tPlayer::inherit(tPlayer *from) {
   ancestor=from;
   ancestor->nrPointingAtMe++;
   if (g::deterministic) {
      for(int i=0;i<genes;i++){
         if(ppc[i] < 0.0) {
            probs[i]=from->inheritGene(i); // inherits w/ chance of mutation
         } else {
            probs[i]=ppc[i];
         }
      }
		for(int i=0;i<genes;i++)
			if(randDouble()<probs[i])
				action=(action<<1)+1;
			else
				action=action<<1;
   } else { /// probabilistic
      for(int i=0;i<genes;i++){
         if(ppc[i] < 0.0) {
            probs[i]=from->inheritGene(i); // inherits w/ chance of mutation
         } else {
            probs[i]=ppc[i];
         }
      }
   }
}

double tPlayer::inheritGene(int w) {
   if(randDouble()<g::mutationRate) {
		if (g::deterministic) return floor(randDouble()+0.5);
      else return randDouble();
	} else
      return probs[w];
}

void tPlayer::exportLOD(FILE *f) {
   if(ancestor!=NULL)
      ancestor->exportLOD(f);
   if (g::savePlays)
      fprintf( f, "%s %d %d %f %f %f %f %f %f %f %f %f %f %s\n",
            cstr(g::experimentID), g::replicateID, born,
            g::rMultiplier, g::zeta, g::beta, g::gamma,
            probs[0], probs[1],
            (1.0-probs[0])*probs[1],
            (1.0-probs[0])*(1.0-probs[1]),
            probs[0]*probs[1],
            probs[0]*(1.0-probs[1]),
            getPlayHistoryAsString()
				);
   else
      fprintf( f, "%s %d %d %f %f %f %f %f %f %f %f %f %f\n",
            cstr(g::experimentID), g::replicateID, born,
            g::rMultiplier, g::zeta, g::beta, g::gamma,
            probs[0], probs[1],
            (1.0-probs[0])*probs[1],
            (1.0-probs[0])*(1.0-probs[1]),
            probs[0]*probs[1],
            probs[0]*(1.0-probs[1])
			  	);
}

char* tPlayer::getPlayHistoryAsString(){
   char* history = new char[plays.size()+1];
   history[plays.size()] = '\0';
   for (int play_i=plays.size()-1; play_i>=0; play_i--) {
      history[play_i] = playToChar[(int)plays[play_i]];
   }
   return history;
}
	
