//
//  main.cpp
//  PTG
//
//  Created by Arend Hintze on 4/11/14.
//  Copyright (c) 2014 Arend Hintze. All rights reserved.
//  Modified by Jory Schossau 2015
//
#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <map>
#include <cmath> // for pow
#include <ctime> // for time
#include <iostream> // for console display
#include <thread> // for multithreading
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

double perSiteMutationRate=0.005;
//int g::update=0;
size_t repeats=1;
int maxAgent=1024;
int totalGenerations=200;
//const int mp=1;

enum class Selection : int {
		task = 1,
		phi = 2,
		r = 4,
		topology = 8,
		genome = 16,
		none = 32
};

map<string, Selection> SelectionValues{
 make_pair("task",Selection::task),
 make_pair("phi",Selection::phi),
 make_pair("r",Selection::r),
 make_pair("topology",Selection::topology),
 make_pair("genome",Selection::genome),
 make_pair("none",Selection::none)
};

bool isaSelectionValue(string value) {
	if (SelectionValues.find(value) == SelectionValues.end())
		return false;
	else
		return true;
}

bool fileExists(const string& name) {
   struct stat buffer;
   return( stat( name.c_str(), &buffer) == 0); // 0 if no problem and file is existing
}

//void computeLOD(FILE *f,FILE *g, tAgent *agent,tGame *game);
const char* cstr(string s) { return s.c_str(); }
double randDouble() { return ((double)rand() / (double)RAND_MAX); }

// processes a subsection of the population
/*void threadedEvaluateFitness(int chunkBegin, int chunkEnd, const vector<tAgent*>& agents, tGame& game, int& evaluations) {
	for (int chunk_index=chunkBegin; chunk_index<chunkEnd; chunk_index++) {
		for (int replicate_i=0; replicate_i < evaluations; replicate_i++) {
			game.executeGame(agents[chunk_index], 2, nullptr, true, -1, -1);
			agents[chunk_index]->fitnesses.push_back(agents[chunk_index]->fitness);
		}
	}
}*/

#ifdef NEVER
/*
int oldmain(int argc, char *argv[]) {
	string g::experimentID;
	int g::replicateID=0;
	vector<tAgent*>agent;
	vector<tAgent*>nextGen;
	int who=0;
	int i;
	tAgent *masterAgent=nullptr;
	tGame *game=nullptr;
	FILE *LODFile=nullptr;
	FILE *genomeFile=nullptr;
	bool showhelp;
	//float evolvePhiLimit;
	//int evolvePhiGenLimit;
	//float evolveRLimit;
	//int evolveRGenLimit;
	//float evolveTopologyLimit;
	//int evolveTopologyGenLimit;
	//int evolveGenomeSizeLimit;
	//int evolveGenomeSizeGenLimit;
	float regimeValueLimit; // evolve until value of fitness is reached, then drop into basic task fitness
	float maxFitness; // holds max fitness of population, used for proportional selection
	float tempFitness; // holds max/actual fitness while being geometrically accumulated
	int regimeGenLimit; // evolve until evolveGenLimit generations passed, then drop into basic task fitness
	vector<string> selectionOn;
	int startGenes;
	bool stopOnLimit;
	int selectionRegime = 0; // bitmask for how to perform selection (fitness, phi, r, topology, genome)
	int nregimes = 0;
	int evaluations=0;						// how many times to test evaluating fitness of an agent
	string filenameLOD, filenameGenome, filenameStartWith;
	int nthreads=0;
	vector<thread> threads;
   //int nknockouts=0; /// defined globally above
   /// Binomial Distribution for determining which locations are mutated

	addp(TYPE::BOOL, &showhelp);
	addp(TYPE::STRING, &filenameLOD, "--lod", "filename to save Line of Descent.");
	addp(TYPE::STRING, &g::experimentID, "--experiment", "unique identifier for this experiment, shared by all replicates.");
	addp(TYPE::INT, &g::replicateID, "--replicate", "unique number to identify this replicate in this experiment.");
	addp(TYPE::STRING, &filenameGenome, "--genome", "filename to save genomes of the LODFile.");
	addp(TYPE::INT, &totalGenerations, "200", false, "--generations", "number of generations to simulate (g::updates).");
	addp(TYPE::STRING, &filenameStartWith, "none", false, "--startWith", "specify a genome file used to seed the population.");
	addp(TYPE::BOOL, &stopOnLimit, "false", false, "--stopOnLimit", "if a limit is specified, then the simulation will stop at the limit.");
	addp(TYPE::INT, &startGenes, "5000", false, "--startGenes", "number of genes with which first organisms begin.");
	addp(TYPE::INT, &evaluations, "1", false, "--evaluations", "number of evaluations for an agent's fitness evaluation.");
	addp(TYPE::INT, &nthreads, "1", false, "--nthreads", (string("number of threads to use. This system reports ")+to_string(thread::hardware_concurrency())+string(" cores available.")).c_str());
	addp(TYPE::STRING, &selectionOn, -1, false, "--selectionOn", "list of parameters on which to perform selection. Valid options are: task, phi, r, topology, genome, none.");
	addp(TYPE::INT, &regimeGenLimit, "-1", false, "--regimeGenLimit", "generation limit to use in all regimes before going back to task fitness.");
	addp(TYPE::FLOAT, &regimeValueLimit, "-1", false, "--regimeValueLimit", "value limit to use in all regimes before going back to task fitness.");
   addp(TYPE::INT, &nknockouts, "0", false, "--nknockouts", "number of knockouts to perform (all n-way knockouts).");
   addp(TYPE::INT, &nstateknockouts, "0", false, "--nstateknockouts", "number of state knockouts to perform (all n-way knockouts).");
   addp(TYPE::BOOL, &phiKnockouts, "false", false, "--phiKnockouts", "Whether to inundate system with noise during analysis (sweeps 0.0-1.0");

	argparse(argv);
	if (showhelp) {
		std::cout << argdetails() << std::endl;
		std::cout << "Example minimal invocation:" << std::endl;
		std::cout << argv[0] << " --experiment=linearFitness --replicate=1 --lod=lineOfDescent.lod --genome=genome.gen" << std::endl;
		std::cout << "or" << std::endl;
		std::cout << argv[0] << " --experiment linearFitness --replicate 1 --lod lineOfDescent.lod --genome genome.gen" << std::endl;
		std::cout << std::endl;
		exit(0);
	}

	/// inputs for selectionRegime from --selectionOn
	if (selectionOn.size() == 0) {
		selectionOn.push_back("task");
		//std::cout << "no regimes specified for selection after --selectionOn option." << std::endl;
		//std::cout << std::endl;
		//exit(1);
	}
	for (string& argument : selectionOn) {
		if (isaSelectionValue(argument)) {
			selectionRegime |= (int)SelectionValues[argument];
		} else {
			std::cout << "invalid --selectionOn option: '" << argument << "'" << std::endl;
			std::cout << std::endl;
			exit(1);
		}
	}

	/// count number of (unique & valid) regimes specified
	if ((selectionRegime & (int)Selection::none) == (int)Selection::none) {
		nregimes = 1;
		selectionRegime = (int)Selection::none;
	} else {
		nregimes = 0;
		int regimeBitmask = selectionRegime;
		for (nregimes = 0; regimeBitmask; nregimes++) { 
			regimeBitmask &= regimeBitmask - 1;
		}
	}

    srand(getpid());
    masterAgent=new tAgent();

	LODFile=fopen(cstr(filenameLOD),"w+t");
	genomeFile=fopen(cstr(filenameGenome),"w+t");	
	srand(getpid());
	agent.resize(maxAgent);
	masterAgent=new tAgent;
	vector<vector<int> > data;
	game=new tGame;
	masterAgent->setupRandomAgent(startGenes);
	
	masterAgent->setupPhenotype();

   float loadingMutationRate = 0.0f;
	if (filenameStartWith != "none") {
		masterAgent->loadAgent(filenameStartWith.c_str());
      loadingMutationRate = 0.0;
	} else {
      loadingMutationRate = 0.5;
	}
   for (i=maxNodes-1; i>=0; --i)
      nodeMask[i]=false;
   if (phiKnockouts) {
      std::cout << "performing phi knockout" << std::endl;
      FILE* koresults = nullptr;
      string koresults_name = g::experimentID+".ko";

      // error checking
      if (filenameStartWith == "none") {
         printf("Error: must provide genome file to test for state knockouts. (--startWith somefile.gen)\n");
         return(1);
      }
      // ok, error checking out of the way. Now for the actual code...
      masterAgent->loadAgent(filenameStartWith.c_str());
      agent[0]=masterAgent;

      // write data
      if (fileExists(koresults_name) == false) {
         koresults = fopen(cstr(koresults_name), "a+t");
         fprintf(koresults, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "g::experimentID", "g::replicateID", "fitness", "phi", "r", "topology", "correct", "incorrect", "knockouts", "nodeids", "noise_level");
      } else {
         koresults = fopen(cstr(koresults_name), "a+t");
      }
      int ko_i=0;
      /// Do noise application here <<<
      for (int eval_i=0; eval_i<evaluations; ++eval_i) {
         float phiKnockoutLevel = 0.0f;
         // 0%-level knockout
         for (ko_i=maxNodes-1; ko_i>=0; --ko_i) {
            nodeMask[ko_i]=false;
         }
         agent[0]->correct=0;
         agent[0]->incorrect=0;
         agent[0]->Phi=0.0;
         agent[0]->R=0.0;
         agent[0]->Topology=0.0;
         agent[0]->NewBinomialDistribution(maxNodes, phiKnockoutLevel);
         threadedEvaluateFitness(0, 1, ref(agent), ref(*game), ref(evaluations));
         tempFitness = 1.0f;
         tempFitness *= pow(fitnessPower,agent[0]->correct - agent[0]->incorrect);
         tempFitness += 1.0f;
         agent[0]->fitness = tempFitness;
         fprintf(koresults, 
               "%s\t%i\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%s\t%f\n",
               cstr(g::experimentID), 
               g::replicateID,
               agent[0]->fitness, 
               agent[0]->Phi, 
               agent[0]->R, 
               agent[0]->Topology, 
               agent[0]->correct, 
               agent[0]->incorrect, 
               0, 
               "-1", 
               phiKnockoutLevel);
         // 10%-level knockout
         phiKnockoutLevel = 0.1f;
         agent[0]->correct=0;
         agent[0]->incorrect=0;
         agent[0]->Phi=0.0;
         agent[0]->R=0.0;
         agent[0]->Topology=0.0;
         agent[0]->NewBinomialDistribution(maxNodes, phiKnockoutLevel);
         threadedEvaluateFitness(0, 1, ref(agent), ref(*game), ref(evaluations));
         tempFitness = 1.0f;
         tempFitness *= pow(fitnessPower,agent[0]->correct - agent[0]->incorrect);
         tempFitness += 1.0f;
         agent[0]->fitness = tempFitness;
         fprintf(koresults, 
               "%s\t%i\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%s\t%f\n",
               cstr(g::experimentID), 
               g::replicateID,
               agent[0]->fitness, 
               agent[0]->Phi, 
               agent[0]->R, 
               agent[0]->Topology, 
               agent[0]->correct, 
               agent[0]->incorrect, 
               0, 
               "-1", 
               phiKnockoutLevel);
         // 1%-9% knockouts
         for (int knockout_i=1; knockout_i<10; ++knockout_i) { // do 1-9 inclusive
            phiKnockoutLevel = (float)knockout_i / (float)100;
            agent[0]->correct=0;
            agent[0]->incorrect=0;
            agent[0]->Phi=0.0;
            agent[0]->R=0.0;
            agent[0]->Topology=0.0;
            agent[0]->NewBinomialDistribution(maxNodes, phiKnockoutLevel);
            threadedEvaluateFitness(0, 1, ref(agent), ref(*game), ref(evaluations));
            tempFitness = 1.0f;
            tempFitness *= pow(fitnessPower,agent[0]->correct - agent[0]->incorrect);
            tempFitness += 1.0f;
            agent[0]->fitness = tempFitness;
            fprintf(koresults, 
                  "%s\t%i\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%s\t%f\n", 
                  cstr(g::experimentID), 
                  g::replicateID, 
                  agent[0]->fitness, 
                  agent[0]->Phi, 
                  agent[0]->R, 
                  agent[0]->Topology, 
                  agent[0]->correct, 
                  agent[0]->incorrect, 
                  0, 
                  "-1", 
                  phiKnockoutLevel);
         }
      }
      fclose(koresults);
      return(0);
   } 
   if (nstateknockouts != 0) {
      FILE* koresults = nullptr;
      string koresults_name = g::experimentID+".ko";

      // error checking
      if (filenameStartWith == "none") {
         printf("Error: must provide genome file to test for state knockouts. (--startWith somefile.gen)\n");
         return(1);
      }
      if (nstateknockouts < 0) {
         printf("Error: --nstateknockouts must be non-negative.\n");
      }

      // ok, error checking out of the way. Now for the actual code...
      masterAgent->loadAgent(filenameStartWith.c_str());
      if (maxNodes < nknockouts) {
         fprintf(stderr, "Agent contains not enough states (defined by maxNodes in globalConst).\n");
         return(1);
      }
      agent[0]=masterAgent;

      // write data
      if (fileExists(koresults_name) == false) {
         koresults = fopen(cstr(koresults_name), "a+t");
         fprintf(koresults, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "g::experimentID", "g::replicateID", "fitness", "phi", "r", "topology", "correct", "incorrect", "knockouts", "nodeids", "noise_nevel");
      } else {
         koresults = fopen(cstr(koresults_name), "a+t");
      }
      bool* mask = new bool[maxNodes];
      int ko_i=0;
      int posOfNextOrder=0;
      bool movedAbit=true;
      for (ko_i=0; ko_i<maxNodes; ++ko_i) {
         mask[ko_i] = false;
      }
      for (ko_i=maxNodes-1; ko_i >= (int)(maxNodes-nstateknockouts); --ko_i) {
         mask[ko_i] = true; // initialize  bits
      }
      // do initial evaluation
      threadedEvaluateFitness(0, 1, ref(agent), ref(*game), ref(evaluations));
      tempFitness = 1.0f;
      tempFitness *= pow(fitnessPower,agent[0]->correct - agent[0]->incorrect);
      tempFitness += 1.0f;
      agent[0]->fitness = tempFitness;
      fprintf(koresults, "%s\t%i\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%s\t%f\n", cstr(g::experimentID), g::replicateID, agent[0]->fitness, agent[0]->Phi, agent[0]->R, agent[0]->Topology, agent[0]->correct, agent[0]->incorrect, 0, "-1", -1.0f);

      {
         // do knockout evaluations
         do {
            string states_string("");
            for (ko_i=maxNodes-1; ko_i>=0; --ko_i) {
               nodeMask[ko_i]=mask[ko_i];
               if (mask[ko_i]) { // build string record of which nodes are toggled off
                  states_string+=std::to_string(ko_i)+",";
               }
            }
            states_string = states_string.substr(0, states_string.size()-1); // remove last ','
            agent[0]->correct=0;
            agent[0]->incorrect=0;
            agent[0]->Phi=0.0;
            agent[0]->R=0.0;
            agent[0]->Topology=0.0;
            threadedEvaluateFitness(0, 1, ref(agent), ref(*game), ref(evaluations));
            tempFitness = 1.0f;
            tempFitness *= pow(fitnessPower,agent[0]->correct - agent[0]->incorrect);
            tempFitness += 1.0f;
            agent[0]->fitness = tempFitness;
            fprintf(koresults, "%s\t%i\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%s\t%f\n", cstr(g::experimentID), g::replicateID, agent[0]->fitness, agent[0]->Phi, agent[0]->R, agent[0]->Topology, agent[0]->correct, agent[0]->incorrect, nstateknockouts, states_string.c_str(), 0.0f);

            // go to next configuration
            movedAbit=false; // start with assumption: all bits are already in lowest order so no bits can be moved to lower order
            for (ko_i=1; ko_i<maxNodes; ++ko_i) {
               if (mask[ko_i] && !mask[ko_i-1]) { // if bit has an empty neighbor (not in lowest order), then we can move the bit
                  mask[ko_i]=false;
                  mask[ko_i-1]=true;
                  movedAbit=true;
                  if ((ko_i > 1) && (mask[0])) { // if lowest bit was already in lowest order
                     posOfNextOrder=ko_i-1; // where this bit was just moved to
                     if (mask[posOfNextOrder-1] == false) { // if next lower bit is not already in its new reset location
                        for (ko_i=posOfNextOrder-1; ko_i>=0; --ko_i) {
                           if (mask[ko_i]) {
                              mask[ko_i]=false;
                              mask[--posOfNextOrder]=true;
                           }
                        }
                     }
                  }
                  break;
               }
            }
         } while (movedAbit);
      }
      delete[] mask;
      fclose(koresults);
      //printf("from main: %p\n", (void *)nodeMask);
      return(0);
   }
   if (nknockouts != 0) {
      FILE* koresults = nullptr;
      string koresults_name = g::experimentID+".ko";

      // error checking
      if (filenameStartWith == "none") {
         printf("Error: must provide genome file to test for knockouts. (--startWith somefile.gen)\n");
         return(1);
      }
      if (nknockouts < 0) {
         printf("Error: --nknockouts must be non-negative.\n");
      }

      // ok, error checking out of the way. Now for the actual code...
      masterAgent->loadAgent(filenameStartWith.c_str());
      int ngates = masterAgent->countGates();
      if (ngates == 0) {
         fprintf(stderr, "Agent contains no gates.\n");
         return(1);
      }
      if (ngates < nknockouts) {
         fprintf(stderr, "Agent contains not enough gates.\n");
         return(1);
      }
      agent[0]=masterAgent;

      // write data
      if (fileExists(koresults_name) == false) {
         koresults = fopen(cstr(koresults_name), "a+t");
         fprintf(koresults, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "g::experimentID", "g::replicateID", "fitness", "phi", "r", "topology", "correct", "incorrect", "knockouts", "nodeids", "noiselevel");
      } else {
         koresults = fopen(cstr(koresults_name), "a+t");
      }
      bool* mask = new bool[agent[0]->hmmus.size()];
      int ko_i=0;
      int posOfNextOrder=0;
      bool movedAbit=true;
      for (ko_i=0; ko_i<(int)agent[0]->hmmus.size(); ++ko_i) {
         mask[ko_i] = false;
      }
      for (ko_i=agent[0]->hmmus.size()-1; ko_i >= (int)agent[0]->hmmus.size()-nknockouts; --ko_i) {
         mask[ko_i] = true; // initialize  bits
      }
      // do initial evaluation
      threadedEvaluateFitness(0, 1, ref(agent), ref(*game), ref(evaluations));
      tempFitness = 1.0f;
      tempFitness *= pow(fitnessPower,agent[0]->correct - agent[0]->incorrect);
      tempFitness += 1.0f;
      agent[0]->fitness = tempFitness;
      fprintf(koresults, "%s\t%i\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%s\t%f\n", cstr(g::experimentID), g::replicateID, agent[0]->fitness, agent[0]->Phi, agent[0]->R, agent[0]->Topology, agent[0]->correct, agent[0]->incorrect, 0, "-1", -1.0f);

      // do knockout evaluations
      do {
         
         string states_string("");
         for (ko_i=agent[0]->hmmus.size()-1; ko_i>=0; --ko_i) {
            agent[0]->hmmus[ko_i]->knockedOut = mask[ko_i];
            if (mask[ko_i]) { // build string record of which nodes are toggled off
               states_string+=std::to_string(ko_i)+",";
            }
         }
         states_string = states_string.substr(0, states_string.size()-1); // remove last ','
         agent[0]->correct=0;
         agent[0]->incorrect=0;
         agent[0]->Phi=0.0;
         agent[0]->R=0.0;
         agent[0]->Topology=0.0;
         threadedEvaluateFitness(0, 1, ref(agent), ref(*game), ref(evaluations));
         tempFitness = 1.0f;
         tempFitness *= pow(fitnessPower,agent[0]->correct - agent[0]->incorrect);
         tempFitness += 1.0f;
         agent[0]->fitness = tempFitness;
         fprintf(koresults, "%s\t%i\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%s\t%f\n", cstr(g::experimentID), g::replicateID, agent[0]->fitness, agent[0]->Phi, agent[0]->R, agent[0]->Topology, agent[0]->correct, agent[0]->incorrect, nknockouts, states_string.c_str(), -1.0f);

         // go to next configuration
         movedAbit=false; // start with assumption: all bits are already in lowest order so no bits can be moved to lower order
         if (ngates != 1) {
            for (ko_i=1; ko_i<(int)agent[0]->hmmus.size(); ++ko_i) {
               if (mask[ko_i] && !mask[ko_i-1]) { // if bit has an empty neighbor (not in lowest order), then we can move the bit
                  mask[ko_i]=false;
                  mask[ko_i-1]=true;
                  movedAbit=true;
                  if ((ko_i > 1) && (mask[0])) { // if lowest bit was already in lowest order
                     posOfNextOrder=ko_i-1; // where this bit was just moved to
                     if (mask[posOfNextOrder-1] == false) { // if next lower bit is not already in its new reset location
                        for (ko_i=posOfNextOrder-1; ko_i>=0; --ko_i) {
                           if (mask[ko_i]) {
                              mask[ko_i]=false;
                              mask[--posOfNextOrder]=true;
                           }
                        }
                     }
                  }
                  break;
               }
            }
         }
      } while (movedAbit);
      delete[] mask;
      fclose(koresults);
      return(0);
   }
   // create population
   for(i=0;i<(int)agent.size();i++){
      agent[i]=new tAgent;
      agent[i]->inherit(masterAgent,loadingMutationRate,0); // large mutation rate to spread the genotypic population
   }
	nextGen.resize(agent.size());
	masterAgent->nrPointingAtMe--;
	std::cout<<"setup complete"<<std::endl;
	printf("%s	%s	%s	%s	%s	%s	%s\n", "g::update","(double)maxFitness","maxPhi","r", "maxTopology", "agent[who]->correct","agent[who]->incorrect");

	while(g::update<totalGenerations) {
		for(i=0;i<(int)agent.size();i++) {
			agent[i]->fitness=0.0;
			agent[i]->fitnesses.clear();
		}

		/// perform fitness evaluation if not none selection specified
		if ((selectionRegime & (int)Selection::none) != (int)Selection::none) {
			threads.clear();
			int chunksize=agent.size()/nthreads;
			for (int threadid=0; threadid < nthreads; threadid++)
				threads.push_back(thread(threadedEvaluateFitness, chunksize*threadid, chunksize*threadid+chunksize, ref(agent), ref(*game), ref(evaluations)));
			if (agent.size()%nthreads != 0) // take care of any uneven division of workload
			{
				threads.push_back(thread(threadedEvaluateFitness, nthreads*chunksize, agent.size(), ref(agent), ref(*game), ref(evaluations)));
			}
			for (thread& t : threads) t.join(); // wait for all threads to finish
		}

		/// perform selection
		maxFitness = 0.0f;
		for(i=0;i<agent.size();i++) {
			tempFitness = 1.0f;
			if ((selectionRegime & (int)Selection::task) == (int)Selection::task) {
				tempFitness *= pow(fitnessPower,agent[i]->correct - agent[i]->incorrect);
				tempFitness += 1.0f;
			}
			if ((selectionRegime & (int)Selection::phi) == (int)Selection::phi) {
				tempFitness *= agent[i]->Phi/16.0f + 1.0f;
			}
			if ((selectionRegime & (int)Selection::r) == (int)Selection::r) {
				tempFitness *= agent[i]->R/4.0f + 1.0f;
			}
			if ((selectionRegime & (int)Selection::topology) == (int)Selection::topology) {
				tempFitness *= (agent[i]->Topology / 100.0f) + 1.0f;
			}
			if ((selectionRegime & (int)Selection::genome) == (int)Selection::genome) {
				tempFitness *= (agent[i]->genome.size() / 20000.0f) + 1.0f;
			}
			if ((selectionRegime & (int)Selection::none) == (int)Selection::none) {
				tempFitness = 2.0f;
				agent[i]->correct = 0;
				agent[i]->incorrect = 0;
				agent[i]->fitness = 1.0f;
				agent[i]->R = 1.0f;
				agent[i]->Phi = 1.0f;
			}
			agent[i]->fitness = tempFitness;
			if (tempFitness > maxFitness) {
				maxFitness = tempFitness;
				who = i;
			}
		}
		if ((selectionRegime & (int)Selection::none) == (int)Selection::none) {
			maxFitness = -1.0f;
		}

		/// regime switch conditions
		if ((nregimes > 1) || (selectionRegime != (int)Selection::task)) {
			if ((regimeGenLimit != -1) && (g::update > regimeGenLimit)) {
				selectionRegime = (int)Selection::task;
				nregimes = 1;
			}
			if (regimeValueLimit > 0.0f) {
				if (maxFitness > pow(fitnessPower, regimeValueLimit*80.0f)) {
					selectionRegime = (int)Selection::task;
					nregimes = 1;
					if (stopOnLimit) break;
				}
			}
		}
			// convert maxValue to 0-100
		printf("%i	%f	%f	%f	%f	%i	%i\n", g::update, maxFitness, agent[who]->Phi, agent[who]->R, agent[who]->Topology, agent[who]->correct, agent[who]->incorrect);

		int j=0;
		for(i=0;i<agent.size();i++) {
			tAgent *d;
			d=new tAgent;
			if(maxFitness<0.0f){
				j=rand()%(int)agent.size();
			} else {
				do{
					j=rand()%(int)agent.size();
				} while((j==(i))||(randDouble>( agent[j]->fitness / maxFitness )));
			}
			d->inherit(agent[j],perSiteMutationRate,g::update);
			nextGen[i]=d;
		}
	//}
		for(i=0;i<agent.size();i++){
			agent[i]->retire();
			agent[i]->nrPointingAtMe--;
			if(agent[i]->nrPointingAtMe==0)
				delete agent[i];
			agent[i]=nextGen[i];
		}
		agent=nextGen;
		g::update++;
	}
*/
#endif

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
//float expx[512];
//double pcd[2]={-1.0,-1.0};
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

struct ScoreDelta {
	unsigned int x,y;
	double value;
	char play;
};

// processes a subsection of the population
void threadedEvaluateFitness(int chunkBegin, int chunkEnd, const vector<vector<tPlayer*>>& player, PopulationFrequencies& frequencies, vector<ScoreDelta>& scoreDeltas) {
	int row(0), col(0), z(0);
   double pool(0.0);
	unsigned char N[neighbors], done[neighbors+1];
	frequencies.c=0; frequencies.d=0; frequencies.m=0; frequencies.i=0; 
	for (int chunk_index=chunkBegin; chunk_index<chunkEnd; chunk_index++) {
		row = chunk_index / xDim;
		col = chunk_index % xDim;
		for(z=0;z<neighbors;z++) {
			N[z] = 0;
		}
		for(z=0;z<neighbors+1;z++){ // for all neighbors and self as last 'neighbor'
			done[z]=player[(row+xm[z])&(xDim-1)][(col+ym[z])&(yDim-1)]->move();
			//player[(row+xm[z])&(xDim-1)][(col+ym[z])&(yDim-1)]->plays.push_back(done[z]); // record play
			N[done[z]]++;
		}
      frequencies.c+=N[C]; frequencies.d+=N[D]; frequencies.m+=N[M]; frequencies.i+=N[I];
		//if ( randDouble() <  expx[ N[C]+N[M]+1 - g::rMultiplier + neighbors+1 ] ) pool = g::rMultiplier;
		//if ( randDouble() < 1.0 / ( 1.0 + exp(-( N[C]+N[M]+1 - g::rMultiplier )) ) ) pool = g::rMultiplier;
		if ( randDouble() < 1.0 / ( 1.0 + exp(-( N[C]+N[M]+1 - g::rMultiplier )) ) ) pool = neighbors+1;
		else pool = 0.0;
		for(z=0;z<neighbors+1;z++) {
			ScoreDelta score;
			score.play = done[z];
			score.x=(row+xm[z])&(xDim-1);
			score.y=(col+ym[z])&(yDim-1);
			switch(done[z]){
				case C: //C ooperator
					//player[(row+xm[z])&(xDim-1)][(col+ym[z])&(yDim-1)]->score+=( ((g::rMultiplier*pool)/((double)neighbors+1.0))-1.0);//offset_min )/offset_max;
					score.value=( ((g::rMultiplier*pool)/((double)neighbors+1.0))-1.0);
					break;
				case D: //D efector
					//player[(row+xm[z])&(xDim-1)][(col+ym[z])&(yDim-1)]->score+=( ((g::rMultiplier*pool)/((double)neighbors+1.0))-(g::beta*((double)N[M]+(double)N[I])/((double)neighbors)));//offset_min )/offset_max;
					score.value=( ((g::rMultiplier*pool)/((double)neighbors+1.0))-(g::beta*((double)N[M]+(double)N[I])/((double)neighbors)));
					break;
				case M://M oralist
					//player[(row+xm[z])&(xDim-1)][(col+ym[z])&(yDim-1)]->score+=( ((g::rMultiplier*pool)/((double)neighbors+1.0))-1.0-(g::gamma*((double)N[D]+(double)N[I])/((double)neighbors)));//offset_min )/offset_max;
					score.value=( ((g::rMultiplier*pool)/((double)neighbors+1.0))-1.0-(g::gamma*((double)N[D]+(double)N[I])/((double)neighbors)));
					break;
				case I://I moralist
					//player[(row+xm[z])&(xDim-1)][(col+ym[z])&(yDim-1)]->score+=( ((g::rMultiplier*pool)/((double)neighbors+1.0))-(g::beta*((double)N[M]+(double)N[I]-(double)1.0)/((double)neighbors))-(g::gamma*((double)N[D]+(double)N[I])/((double)neighbors)));//offset_min )/offset_max;
					score.value=( ((g::rMultiplier*pool)/((double)neighbors+1.0))-(g::beta*((double)N[M]+(double)N[I]-(double)1.0)/((double)neighbors))-(g::gamma*((double)N[D]+(double)N[I])/((double)neighbors)));
					break;
			}
			//player[(row+xm[z])&(xDim-1)][(col+ym[z])&(yDim-1)]->score = fmax(0.0,player[(row+xm[z])&(xDim-1)][(col+ym[z])&(yDim-1)]->score);
			//scoreDeltas.push_back(score);
			scoreDeltas[chunk_index-chunkBegin] = score;
		}
	}
}

int main (int argc, char* argv[]) {
   int x,y,z,i;
   double maxFit,localFit;
   unsigned int maxFitPlayerX,maxFitPlayerY;
   vector< vector<tPlayer*> > player;
   vector< PopulationFrequencies* > census;
	PopulationFrequencies* frequencies;
   int lxm[neighbors],lym[neighbors];
   srand(getpid());
   bool debug;

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
   int nthreads;
	vector<thread> threads;

   
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
	addp(TYPE::INT, &nthreads, "1", false, "--nthreads", (string("number of threads to use. This system reports ")+ to_string(thread::hardware_concurrency())+ string(" cores available.")).c_str());
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
		std::cout << std::endl;
		std::cout << string("There are ") + to_string(neighbors) + string(" neighbors defined in this build.") << std::endl;
		std::cout << std::endl;
		exit(0);
	}
   float bounty = 0.0;
   float offset_min = 0.0;
   float offset_max = 0.0;
	// precompute the exponential function of # coop and probability of success
//	for (int i= 0; i<2*(neighbors+1); i++) { /// remember to offset by neighbors when indexing into expx ex: expx[2+neighbors]
//		expx[i] = 1.0 / ( 1.0 + exp(-i+neighbors+1) );
//	}
//   //C ooperator
//   ((g::rMultiplier*pool)/((double)neighbors+1.0))-1.0;
//   //D efector
//   ((g::rMultiplier*pool)/((double)neighbors+1.0))-(g::beta*((double)N[M]+(double)N[I])/((double)neighbors));
//   //M oralist
//   ((g::rMultiplier*pool)/((double)neighbors+1.0))-1.0-(g::gamma*((double)N[D]+(double)N[I])/((double)neighbors));
//   //I moralist
//   ((g::rMultiplier*pool)/((double)neighbors+1.0))-(g::beta*((double)N[M]+(double)N[I]-(double)1.0)/((double)neighbors))-(g::gamma*((double)N[D]+(double)N[I])/((double)neighbors));
   //
   player.resize(xDim);
   for(x=0;x<xDim;x++) {
      player[x].resize(yDim);
      for(y=0;y<yDim;y++) {
         player[x][y]=new tPlayer;
		}
   }
	vector<PopulationFrequencies> threadFrequencies;
	threadFrequencies.resize(nthreads);
	PopulationFrequencies nonthreadFrequencies;
	vector<vector<ScoreDelta>> threadScoreDeltas;
	threadScoreDeltas.resize(nthreads);
	vector<ScoreDelta> nonthreadScoreDeltas;
   for(g::update=1;g::update<g::updates;g::update++){
      frequencies = new PopulationFrequencies;
      frequencies->c=0; frequencies->d=0; frequencies->m=0; frequencies->i=0; 
      for(x=0;x<neighbors;x++){ // well-mix
         xm[x] = (rand()&31) * (-1*rand()&1);
         ym[x] = (rand()&31) * (-1*rand()&1);
      }
		/// multithreading version
		threads.clear();
		int chunksize=(xDim*yDim)/nthreads;
		/// resize scoreDeltas to be the size of each thread workload
		for (auto& each_scoreset : threadScoreDeltas) {
			each_scoreset.resize(chunksize);
		}
		nonthreadScoreDeltas.resize( (xDim*yDim)%nthreads );
		/// launch threads
		for (int threadid=0; threadid < nthreads; threadid++) {
			threads.push_back( thread( threadedEvaluateFitness, chunksize*threadid, chunksize*threadid+chunksize, ref(player), ref(threadFrequencies[threadid]), ref(threadScoreDeltas[threadid]) ) );
		}
		/// use host thread for any remainder
		if ((xDim*yDim)%nthreads != 0) { // take care of any uneven division of workload
			threads.push_back( thread( threadedEvaluateFitness, nthreads*chunksize, xDim*yDim, ref(player), ref(nonthreadFrequencies), ref(nonthreadScoreDeltas) ) );
		}
		for (thread& t : threads) t.join(); // wait for all threads to finish
		/// play counts
		for (auto& f : threadFrequencies) { /// threaded frequencies
			frequencies->c+=f.c; frequencies->d+=f.d; frequencies->m+=f.m; frequencies->i+=f.i;
		}
		frequencies->c+=nonthreadFrequencies.c; frequencies->d+=nonthreadFrequencies.d; frequencies->m+=nonthreadFrequencies.m; frequencies->i+=nonthreadFrequencies.i; /// nonthreaded frequencies
		frequencies->round();
		/// fitness gains
		for (auto& scoreCollection : threadScoreDeltas) { /// threaded fitness gains
			for (auto& scoring : scoreCollection) {
				player[scoring.x][scoring.y]->score += scoring.value;
				player[scoring.x][scoring.y]->plays.push_back( scoring.play );
			}
		}
		for (auto& scoring : nonthreadScoreDeltas) { /// nonthreaded fitness gains
			player[scoring.x][scoring.y]->score += scoring.value;
			player[scoring.x][scoring.y]->plays.push_back( scoring.play );
		}
		/// fix possible negative fitnesses
		for (x=0; x<xDim; x++)
			for (y=0; y<yDim; y++)
				player[x][y]->score = fmax(0.0, player[x][y]->score);
		/// single-threaded version
      //for(x=0;x<xDim;x++)
      //   for(y=0;y<yDim;y++){
      //      for(z=0;z<possibleMoves;z++)
      //         N[z]=0;
      //      for(z=0;z<neighbors+1;z++){ // for all neighbors and self as last 'neighbor'
      //         done[z]=player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->move();
      //         player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->plays.push_back(done[z]); // record play
      //         N[done[z]]++;
      //      }
      //      frequencies->c+=N[C]; frequencies->d+=N[D]; frequencies->m+=N[M]; frequencies->i+=N[I];
		//		if ( randDouble() <  expx[ N[C]+N[M]+1 - g::rMultiplier ]) pool = neighbors+1;
		//		else pool = 0.0;
      //      for(z=0;z<neighbors+1;z++) {
      //         switch(done[z]){
      //            case C: //C ooperator
      //               player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score+=( ((g::rMultiplier*pool)/((double)neighbors+1.0))-1.0);//offset_min )/offset_max;
      //               break;
      //            case D: //D efector
      //               player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score+=( ((g::rMultiplier*pool)/((double)neighbors+1.0))-(g::beta*((double)N[M]+(double)N[I])/((double)neighbors)));//offset_min )/offset_max;
      //               break;
      //            case M://M oralist
      //               player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score+=( ((g::rMultiplier*pool)/((double)neighbors+1.0))-1.0-(g::gamma*((double)N[D]+(double)N[I])/((double)neighbors)));//offset_min )/offset_max;
      //               break;
      //            case I://I moralist
      //               player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score+=( ((g::rMultiplier*pool)/((double)neighbors+1.0))-(g::beta*((double)N[M]+(double)N[I]-(double)1.0)/((double)neighbors))-(g::gamma*((double)N[D]+(double)N[I])/((double)neighbors)));//offset_min )/offset_max;
      //               break;
      //         }
      //         player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score = fmax(0.0,player[(x+xm[z])&(xDim-1)][(y+ym[z])&(yDim-1)]->score);
      //      }
      //}
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
      //if (debug && (g::update&511) == 511) std::cout<< g::update<<" "<<maxFit<<" "<<player[maxFitPlayerX][maxFitPlayerY]->probs[0]<<" "<<player[maxFitPlayerX][maxFitPlayerY]->probs[1]<<std::endl;
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
         //fprintf( fileLOD,"%s %s %s %s %s %s %s %s %s %s\n","experiment","replicate","born","p0","p1","pc","pd","pm","pi","plays" );
      else
         fprintf( fileLOD, "%s %s %s %s %s %s %s %s %s %s %s %s %s\n",
               "experiment", "replicate", "born",
               "r","zeta","beta","gamma",
               "p0","p1",
               "pc","pd","pm","pi" );
         //fprintf( fileLOD,"%s %s %s %s %s %s %s %s %s\n","experiment","replicate","born","p0","p1","pc","pd","pm","pi" );
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
         fprintf( fileEnd, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
               "experiment", "replicate",
               "r","zeta","beta","gamma",
               "c","d","m","i",
               "p0","p1",
               "pc","pd","pm","pi" );
      }
      tPlayer* LCA = player[0][0]->ancestor->ancestor->ancestor;
      fprintf(fileEnd, "%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
            cstr(g::experimentID), g::replicateID,
            g::rMultiplier, g::zeta, g::beta, g::gamma,
            frequencies->c, frequencies->d, frequencies->m, frequencies->i,
            LCA->probs[0],LCA->probs[1],
            LCA->probs[0]*LCA->probs[1],
            LCA->probs[0]*(1.0-LCA->probs[1]),
            (1.0-LCA->probs[0])*LCA->probs[1],
            (1.0-LCA->probs[0])*(1.0-LCA->probs[1]));
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
   if (g::deterministic) {
		if (fcdmi[0] != -1.0) { // if exact player types specified
			action=0;
			for(i=0;i<genes;i++)
				if(randDouble()<probs[i])
					action=action<<1;
				else
					action=(action<<1)+1;
		} else { // otherwise pick a random action for this player's life
			action=((rand()&1)<<1) + (rand()&1);
		}
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
      /*fprintf( f,"%s %d %d %f %f %f %f %f %f %s\n",cstr(g::experimentID),g::replicateID,born,probs[0],probs[1],
            (probs[0])*(probs[1]),
            (probs[0])*(1.0-probs[1]),
            (1.0-probs[0])*(probs[1]),
            (1.0-probs[0])*(1.0-probs[1]),
            getPlayHistoryAsString() );*/
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
      /*fprintf( f,"%s %d %d %f %f %f %f %f %f\n",cstr(g::experimentID),g::replicateID,born,probs[0],probs[1],
            (probs[0])*(probs[1]),
            (probs[0])*(1.0-probs[1]),
            (1.0-probs[0])*(probs[1]),
            (1.0-probs[0])*(1.0-probs[1]) );*/
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
	
//	agent[0]->ancestor->saveLOD(LODFile,genomeFile, g::experimentID, g::replicateID, -1); // -1 to tell saveLOD to make header for csv
//	if (stopOnLimit) {
//		float maxPhi=0.0;
//		tAgent* bestAgent=nullptr;
//		for (tAgent* a : agent) {
//			if (a->Phi > maxPhi) {
//				maxPhi = a->Phi;
//				bestAgent = a;
//			}
//		}
//		if (bestAgent) {
//			bestAgent->saveGenome(genomeFile);
//		}
//	} else {
//		agent[0]->ancestor->ancestor->saveGenome(genomeFile);
//	}
////	agent[0]->ancestor->saveToDot(argv[3]);
//	agent.clear();
//	nextGen.clear();
//	delete masterAgent;
//	delete game;
//	return 0;
//}

//void computeLOD(FILE *f,FILE *g, tAgent *agent,tGame *game){
	/*vector<vector<int> > table;
	double R,oldR;
	if(agent->ancestor!=NULL)
		computeLOD(f,g,agent->ancestor,game);
	agent->setupPhenotype();
	table=game->executeGame(agent, 2, NULL,false,-1,-1);
	R=game->computeR(table,0);
	oldR=game->computeOldR(table);
	fprintf(f,"%i	%i	%i	%f	%f",agent->ID,agent->correct,agent->incorrect,agent->extra);
	fprintf(f,"\n");
	*/
//}

