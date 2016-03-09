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

#define C 0
#define D 1
#define M 2
#define I 3

#define genes 2
// xDim, yDim assumed powers of 2
#define xDim 32
#define yDim 32
#define possibleMoves 4

// neighbors is assumed power of 2
#define neighbors 19
//int xm[neighbors+1]={0,1,0,-1,0};
//int ym[neighbors+1]={-1,0,1,0,0};
int xm[neighbors+1]={0,1,0,-1, 1,1,-1,-1, 0,2,0,-2, 2,2,-2,-2, 0};
int ym[neighbors+1]={-1,0,1,0, 1,-1,-1,1, -2,0,2,0, 2,-2,-2,2, 0};
vector< float > pcd; // Pc & Pd inheritances: 0.0 <= p <= 1.0 and is constant value otherwise -1 implies usual mutation
vector< float > fcdmi; // fraction of initial C,D,M,I players

namespace g {
	string experimentID;
	int replicateID;

	double mutationRate;
	double replacementRate=0.02;

	double beta;
	double gamma;
	float rMultiplier;
	double zeta;
	int radius=32;
	int updates=42;
	bool savePlays;
	bool deterministic;
	int update=0;
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
   int x,y,z,i;
   double maxFit,localFit;
   unsigned int maxFitPlayerX,maxFitPlayerY;
   vector< vector<tPlayer*> > player;
   vector< PopulationFrequencies* > census; /// record of frequencies over time
	PopulationFrequencies* frequencies;
   int lxm[neighbors],lym[neighbors];
   srand(getpid());
   bool debug;
	double pool(0.0);
	unsigned char N[neighbors], done[neighbors+1];

   bool showhelp;
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
   /// other parameters

   
   /// required
	addp(TYPE::STRING, &g::experimentID, "--experiment", "unique identifier for this experiment, shared by all replicates.");
	addp(TYPE::INT, &g::replicateID, "--replicate", "unique number to identify this replicate in this experiment.");
   /// not required
	addp(TYPE::BOOL, &showhelp);
   addp(TYPE::BOOL, &debug, "false", false, "--debug", "enables verbose printing.");
   addp(TYPE::BOOL, &g::savePlays, "false", false, "--savePlays", "Exports play history in line of descent output.");
   addp(TYPE::BOOL, &g::deterministic, "false", false, "--deterministic", "Forces a deterministic strategy simulation (random initial pop).");
   addp(TYPE::DOUBLE, &g::beta, "0.02", false, "--beta", "Beta parameter: effect of punishment.");
   addp(TYPE::DOUBLE, &g::gamma, "0.02", false, "--gamma", "Gamma parameter: cost of punishment.");
   addp(TYPE::FLOAT, &g::rMultiplier, "15", false, "--r", "R multiplier parameter.");
   addp(TYPE::DOUBLE, &g::zeta, "1.0", false, "--zeta", "Payoff zeta for group success.");
   addp(TYPE::INT, &g::radius, "32", false, "--radius", "Length of one side of square group region, 32 is K=1024 which implies well-mixed since popoulation is 1024.");
	addp(TYPE::STRING, &filenameLOD, "none", false, "--lod", "filename to save Line of Descent.");
	addp(TYPE::STRING, &filenameGenome, "none", false, "--genome", "filename to save LCA genome.");
	addp(TYPE::STRING, &filenamePOP, "none", false, "--pop", "filename to save population frequencies.");
   addp(TYPE::STRING, &filenameEnd, "none", false, "--end", "filename to save best genome and pop frequencies.");
	addp(TYPE::INT, &g::updates, "200", false, "--updates", "number of updates to simulate (not generations).");
   addp(TYPE::FLOAT, &pcd, 2, "-1.0", false, "--pcd", "The probabilities for cooperation and defection (0-1)");
	addp(TYPE::FLOAT, &fcdmi, 4, "-1.0", false, "--fcdmi", "The initial population frequencies for C D M I (0-1) (forces deterministic)");
   addp(TYPE::DOUBLE, &g::mutationRate, "0.05", false, "--mu", "mutation rate (per site)");

	argparse(argv);
	if (showhelp) {
		std::cout << argdetails() << std::endl;
		std::cout << "Example minimal invocation:" << std::endl;
		std::cout << argv[0] << " --experiment=deterministic --replicate=1 --lod=lineOfDescent.lod --genome=genome.gen" << std::endl;
		std::cout << "or" << std::endl;
		std::cout << argv[0] << " --experiment deterministic --replicate 1 --lod lineOfDescent.lod --genome genome.gen" << std::endl;
		std::cout << "or" << std::endl;
		std::cout << argv[0] << " --experiment trial --replicate 1 --end out.end --deterministic --beta 0 --gamma 0 --fcdmi 0.25 0.25 0.25 0.25 --debug --updates 10000 --zeta 0.5 --r 10" << std::endl;
		std::cout << std::endl;
		std::cout << string("There are ") + to_string(neighbors) + string(" neighbors defined in this build.") << std::endl;
		std::cout << std::endl;
		exit(0);
	}
   player.resize(xDim);
   for(x=0;x<xDim;x++) {
      player[x].resize(yDim);
      for(y=0;y<yDim;y++) {
         player[x][y]=new tPlayer;
		}
   }
   for(g::update=1;g::update<g::updates;g::update++){
      frequencies = new PopulationFrequencies;
      frequencies->c=0; frequencies->d=0; frequencies->m=0; frequencies->i=0; 
      for(x=0;x<neighbors;x++){ // well-mix
         xm[x] = (rand()&31) * (-1*rand()&1);
         ym[x] = (rand()&31) * (-1*rand()&1);
      }
		/// fix possible negative fitnesses
		for (x=0; x<xDim; x++)
			for (y=0; y<yDim; y++)
				player[x][y]->score = fmax(0.0, player[x][y]->score);
		/// single-threaded version
      for(x=0;x<xDim;x++)
         for(y=0;y<yDim;y++){
            for(z=0;z<possibleMoves;z++)
               N[z]=0;
            for(z=0;z<neighbors+1;z++){ // for all neighbors and self as last 'neighbor'
               done[z]=player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->move();
               player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->plays.push_back(done[z]); // record play
               N[done[z]]++;
            }
            frequencies->c+=N[C]; frequencies->d+=N[D]; frequencies->m+=N[M]; frequencies->i+=N[I];
				if ( randDouble() < threshold(N[C]+N[M]+1-g::rMultiplier, Threshold::Stepwise) ) pool = neighbors+1;
				else pool = 0.0;
            for(z=0;z<neighbors+1;z++) {
               switch(done[z]){
                  case C: //C ooperator
                     player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score+=( (((g::rMultiplier)*pool)/((double)neighbors+1.0))-1.0);
                     break;
                  case D: //D efector
                     player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score+=( (((g::rMultiplier)*pool)/((double)neighbors+1.0))-(g::beta*((double)N[M]+(double)N[I])/((double)neighbors)));
                     break;
                  case M://M oralist
                     player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score+=( (((g::rMultiplier)*pool)/((double)neighbors+1.0))-1.0-(g::gamma*((double)N[D]+(double)N[I])/((double)neighbors)));
                     break;
                  case I://I moralist
                     player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score+=( (((g::rMultiplier)*pool)/((double)neighbors+1.0))-(g::beta*((double)N[M]+(double)N[I]-(double)1.0)/((double)neighbors))-(g::gamma*((double)N[D]+(double)N[I])/((double)neighbors)));
                     break;
               }
               player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score = fmax(0.0,player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score);
            }
      }
      if (debug && (g::update&511) == 511) {
         maxFit=0.0;
         for(x=0;x<xDim;x++) {
            for(y=0;y<yDim;y++) {
               if(player[x][y]->score>maxFit) {
                  maxFit=player[x][y]->score;
                  maxFitPlayerX=x;
                  maxFitPlayerY=y;
               }
            }
         }
      }
      for(z=0;z<xDim*yDim*g::replacementRate;z++){
         do{
            x=rand()&(xDim-1);
            y=rand()&(yDim-1);
         } while(player[x][y]->born==g::update);
         localFit=0.0;
         
         for(i=0;i<neighbors;i++){
            do{
               lxm[i]=rand()%g::radius;
               lym[i]=rand()%g::radius;
               if((rand()&1)==0) lxm[i]=-lxm[i];
               if((rand()&1)==0) lym[i]=-lym[i];
            } while((lxm[i]==0)&&(lym[i]==0));
         }
         for(i=0;i<neighbors;i++)
            localFit+=player[(x+lxm[i])&(xDim-1)][(y+lym[i])&(yDim-1)]->score;
         if(localFit!=0.0){
            i=rand()%neighbors;
            while( (player[(x+lxm[i])&(xDim-1)][(y+lym[i])&(yDim-1)]->score/localFit)<randDouble() )
               i=rand()%neighbors; // [JDS] I think we could just use an increment through neighbors, instead of rand
         }
         else
            i=rand()&(neighbors-1);
         
         player[x][y]->nrPointingAtMe--;
         if(player[x][y]->nrPointingAtMe==0)
            delete player[x][y];
         player[x][y]=new tPlayer;
         player[x][y]->inherit(player[(x+lxm[i])&(xDim-1)][(y+lym[i])&(yDim-1)]);
      }
      if (debug && (g::update&511) == 511) std::cout<<frequencies->c<<" "<<frequencies->d<<" "<<frequencies->m<<" "<<frequencies->i<<std::endl;
      frequencies->normalize();
      census.push_back(frequencies);
      if (frequencies->c == 1.0 || frequencies->d == 1.0 || frequencies->m == 1.0 || frequencies->i == 1.0) {
         if (debug) printf( "%f %f %f %f\n",frequencies->c, frequencies->d, frequencies->m, frequencies->i );
         break;
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
      fprintf(fileEnd, "%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
            cstr(g::experimentID), g::replicateID,
            g::rMultiplier, g::zeta, g::beta, g::gamma,
            frequencies->c, frequencies->d, frequencies->m, frequencies->i,
            LCA->probs[0],LCA->probs[1],
            LCA->probs[0]*LCA->probs[1],
            LCA->probs[0]*(1.0-LCA->probs[1]),
            (1.0-LCA->probs[0])*LCA->probs[1],
            (1.0-LCA->probs[0])*(1.0-LCA->probs[1]),
				fcdmi[0], fcdmi[1], fcdmi[2], fcdmi[3]);
      fclose(fileEnd);
   }
   return 0;
}

tPlayer::tPlayer() {
   int i;
	if (fcdmi[0] != -1.0) {
		double rnd = randDouble();
		if (rnd <= fcdmi[0]) { probs[0]=1.0; probs[1]=1.0; }
		else if (rnd <= fcdmi[0]+fcdmi[1]) { probs[0]=1.0; probs[1]=0.0; }
		else if (rnd <= fcdmi[0]+fcdmi[1]+fcdmi[2]) { probs[0]=0.0; probs[1]=1.0; }
		else { probs[0]=0.0; probs[1]=0.0; }
	} else {
		for(i=0;i<genes;i++) {
			if (pcd[i] == -1.0) {
				probs[i]=0.5;
			} else {
				probs[i]=pcd[i];
			}
		}
	}
	if (fcdmi[0] != -1.0) { // if exact player types specified
		action=0;
		for(i=0;i<genes;i++)
			if(randDouble()<probs[i])
				action=action<<1;
			else
				action=(action<<1)+1;
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
         r=r<<1;
      else
         r=(r<<1)+1;
   }
   return r;
}
void tPlayer::inherit(tPlayer *from) {
   ancestor=from;
   ancestor->nrPointingAtMe++;
   if (g::deterministic) {
      if (randDouble()<g::mutationRate)
         action=((rand()&1)<<1) + (rand()&1);
      else
         action=from->action;
   } else {
      for(int i=0;i<genes;i++){
         if(pcd[i] < 0.0) {
            probs[i]=from->inheritGene(i);
         } else {
            probs[i]=pcd[i];
         }
      }
   }
}

double tPlayer::inheritGene(int w) {
   if(randDouble()<g::mutationRate && not g::deterministic)
      return randDouble();
   else
      return probs[w];
}

int breakcount=1000;

void tPlayer::exportLOD(FILE *f) {
   if(ancestor!=NULL)
      ancestor->exportLOD(f);
   if (g::savePlays)
      fprintf( f, "%s %d %d %f %f %f %f %f %f %f %f %f %f %s\n",
            cstr(g::experimentID), g::replicateID, born,
            g::rMultiplier, g::zeta, g::beta, g::gamma,
            probs[0],probs[1],
            probs[0]*probs[1],
            probs[0]*(1.0-probs[1]),
            (1.0-probs[0])*probs[1],
            (1.0-probs[0])*(1.0-probs[1]),
            getPlayHistoryAsString() );
   else
      fprintf( f, "%s %d %d %f %f %f %f %f %f %f %f %f %f\n",
            cstr(g::experimentID), g::replicateID, born,
            g::rMultiplier, g::zeta, g::beta, g::gamma,
            probs[0],probs[1],
            probs[0]*probs[1],
            probs[0]*(1.0-probs[1]),
            (1.0-probs[0])*probs[1],
            (1.0-probs[0])*(1.0-probs[1]) );
}

char* tPlayer::getPlayHistoryAsString(){
   char* history = new char[plays.size()+1];
   history[plays.size()] = '\0';
   for (int play_i=plays.size()-1; play_i>=0; play_i--) {
      history[play_i] = playToChar[(int)plays[play_i]];
   }
   return history;
}
	
