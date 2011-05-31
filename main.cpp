#include "PhyloTree.h"
#include "EvolutionModel.h"
#include "Fasta.h"
#include "EASystem.h"
#include "utils.h"
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
            << "[mut. rate] [recomb. rate] [no of sequences] "
            << "[log file path]" << endl;
        return 1;
    }

    unsigned int popSize, elitism, gens, nSequences;
    double mutRate, recRate;
    string  logFilePath;

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
        ss >> nSequences;
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
    logFile << "# random seed: " << seed << endl;
    logFile << "# pop size: " << popSize << endl;
    logFile << "# elitism: " << elitism << endl;
    logFile << "# generations: " << gens << endl;
    logFile << "# mutation rate: " << mutRate << endl;
    logFile << "# recombination rate: " << recRate << endl;
    logFile << "# log file path: " << logFilePath + logFileName << endl;

    // create a random population of trees
    vector<string> randomTrees;
    for (unsigned int i = 0; i < popSize; i++) {
        vector<PhyloTreeNode*> nodes;
        for (unsigned int j = 0; j < nSequences; j++) {
            nodes.push_back(new PhyloTreeNode(convertToString(j), ""));
        }
        PhyloTree t;
        t.buildRandomTree(nodes);
        randomTrees.push_back(PhyloTreeNode::prefixRepresentation(t.getRoot()));
    }

    EASystem<string> testEA(new MutateTree(nSequences, mutRate),
            new RecombineTree(recRate), new RankSelection<string>,
            new PipeFitnessFunc<string>);
    testEA.setLogStream(&logFile);
    testEA.setElitism(elitism);
    testEA.setDebugging(true);
    testEA.setPopulation(randomTrees);

    testEA.runGenerations(gens);

    cerr << "Final generation:\n";
    testEA.exportGenomes(cerr);
    logFile << "\n\nFinal generation (best first):\n";
    testEA.exportGenomes(logFile);

    logFile.close();

    return 0;
}
