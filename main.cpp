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
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " [pop size] [elitism] [gens] "
            << "[mut. rate] [recomb. rate]" << endl;
        return 1;
    }

    unsigned int popSize, elitism, gens;
    double mutRate, recRate;

    // read arguments and assign their values to the EA parameters
    {
        std::stringstream ss;
        for (unsigned int i = 1; i < 6; i++) {
            ss << argv[i];
        }

        ss >> popSize;
        ss >> elitism;
        ss >> gens;
        ss >> mutRate;
        ss >> recRate;
    }

    unsigned int seed = time(NULL);
    srand(seed);

    // create a log file named "log[seed].txt"
    std::ostringstream oss;
    oss << seed;
    std::ofstream logFile;
    string logFileName("log");
    logFileName.append(oss.str()).append(".txt");
    logFile.open (logFileName.c_str());

    // record seed and EA parameters
    logFile << "# random seed: " << seed << endl;
    logFile << "# pop size: " << popSize << endl;
    logFile << "# elitism: " << elitism << endl;
    logFile << "# generations: " << gens << endl;
    logFile << "# mutation rate: " << mutRate << endl;
    logFile << "# recombination rate: " << recRate << endl;

    // create a random population of trees
    vector<string> randomTrees;
    vector<PhyloTreeNode*> nodes;
    for (unsigned int i = 0; i < popSize; i++) {
        nodes = Fasta::readFastaFile("tests/aligned.fasta", 0);
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
