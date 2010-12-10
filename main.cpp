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
using std::string;
using std::vector;

int main(int argc, const char *argv[])
{
    unsigned int seed = time(NULL);

    std::ostringstream oss;
    oss << seed;
    std::ofstream logFile;
    string logFileName("log");
    logFileName.append(oss.str()).append(".txt");
    logFile.open (logFileName.c_str());

    logFile << "# random seed: " << seed << std::endl;
    srand(seed);
    //Kimura* k = new Kimura(10.0);

    vector<string> randomTrees;
    vector<PhyloTreeNode*> nodes;

    unsigned int popSize = 40;
    unsigned int elitism = 4;
    unsigned int gens = 300;
    double mutRate = 0.1;
    double recRate = 0.7;

    logFile << "# pop size: " << popSize << std::endl;
    logFile << "# elitism: " << elitism << std::endl;
    logFile << "# generations: " << gens << std::endl;
    logFile << "# mutation rate: " << mutRate << std::endl;
    logFile << "# recombination rate: " << recRate << std::endl;

    // create a random population of trees
    for (unsigned int i = 0; i < popSize; i++) {
        nodes = Fasta::readFastaFile("tests/aligned.fasta", 0);
        PhyloTree t;
        t.buildRandomTree(nodes);
        randomTrees.push_back(PhyloTreeNode::prefixRepresentation(t.getRoot()));
    }

    EASystem<string> testEA(new MutateTree(nodes.size(), mutRate), new RecombineTree(recRate), new RankSelection<string>, new PipeFitnessFunc<string>);
    testEA.setLogStream(&logFile);
    testEA.setElitism(elitism);
    testEA.setDebugging(true);
    testEA.setPopulation(randomTrees);
    testEA.runGenerations(gens);

    cerr << "Final generation:\n";
    testEA.exportGenomes(cerr);
    testEA.exportGenomes(logFile);

    logFile.close();

    return 0;
}
