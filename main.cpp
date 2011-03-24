#include "PhyloTree.h"
#include "EvolutionModel.h"
#include "Fasta.h"
#include "EASystem.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

#include "EAOperators.h"
#include "TreeOperators.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

int main(int argc, const char *argv[])
{
    if (argc < 8) {
        cerr << "Usage: " << argv[0] << " [pop size] [elitism] [gens] "
            << "[mut. rate] [recomb. rate] [sequence file.fasta] "
            << "[log file path]" << endl;
        return 1;
    }

    unsigned int popSize, elitism, gens;
    double mutRate, recRate;
    string seqFile, logFilePath;

    // read arguments and assign their values to the EA parameters
    {
	    std::stringstream ss;
	    ss << argv[1];
	    ss >> popSize;
	    ss.clear();

	    ss << argv[2];
	    ss >> elitism;
	    ss.clear();

	    ss << argv[3];
	    ss >> gens;
	    ss.clear();

	    ss << argv[4];
	    ss >> mutRate;
	    ss.clear();

	    ss << argv[5];
	    ss >> recRate;
	    ss.clear();

	    ss << argv[6];
	    ss >> seqFile;
	    ss.clear();

	    ss << argv[7];
	    ss >> logFilePath;
	    ss.clear();
    }

    unsigned int seed = time(NULL);
    srand(seed);

    // create a log file named "log[seed].txt"
    std::ostringstream oss;
    oss << seed;
    std::ofstream logFile;
    string logFileName("log");
    logFileName.append(oss.str()).append(".txt");
    logFile.open ((logFilePath+"/"+logFileName).c_str());

    // record seed and EA parameters
    cerr << "# random seed: " << seed << endl;
    cerr << "# pop size: " << popSize << endl;
    cerr << "# elitism: " << elitism << endl;
    cerr << "# generations: " << gens << endl;
    cerr << "# mutation rate: " << mutRate << endl;
    cerr << "# recombination rate: " << recRate << endl;
    cerr << "# log file path: " << logFilePath + logFileName << endl;

    cerr << "aff\n";

    // create a random population of trees
    vector<string> randomTrees;
    vector<PhyloTreeNode*> nodes;
    for (unsigned int i = 0; i < popSize; i++) {
        nodes = Fasta::readFastaFile(seqFile, 0);
        PhyloTree t;
        t.buildRandomTree(nodes);
        randomTrees.push_back(PhyloTreeNode::prefixRepresentation(t.getRoot()));
    }

    EASystem<string> testEA(new MutateTree(nodes.size(), mutRate),
            new RecombineTree(recRate), new RankSelection<string>,
            new PipeFitnessFunc<string>);
    testEA.setLogStream(&logFile);
    testEA.setElitism(elitism);
    testEA.setDebugging(true);
    testEA.setPopulation(randomTrees);

    testEA.runGenerations(gens);

    cerr << "Final generation:\n";
    testEA.exportGenomes(cerr);
    logFile << "Final generation (best first):\n";
    testEA.exportGenomes(logFile);

    logFile.close();

    return 0;
}
