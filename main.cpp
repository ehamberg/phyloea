#include "PhyloTree.h"
#include "EvolutionModel.h"
#include "Fasta.h"
#include "EASystem.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "EAOperators.h"
#include "TreeOperators.h"

using std::cout;
using std::cerr;
using std::string;
using std::vector;

int main(int argc, const char *argv[])
{
    std::ofstream logFile;
    logFile.open ("log.txt");

    unsigned int seed = 1291705571;//time(NULL);

    logFile << "# random seed: " << seed << std::endl;
    srand(seed);
    //Kimura* k = new Kimura(10.0);

    vector<string> randomTrees;
    vector<PhyloTreeNode*> nodes;

    // create a random population of trees
    for (unsigned int i = 0; i < 8; i++) {
        nodes = Fasta::readFastaFile("tests/aligned.fasta", 0);
        PhyloTree t;
        t.buildRandomTree(nodes);
        randomTrees.push_back(PhyloTreeNode::prefixRepresentation(t.getRoot()));
    }

    EASystem<string> testEA(new MutateTree(nodes.size(), 0.1), new RecombineTree(0.7), new RankSelection<string>, new PipeFitnessFunc<string>);
    testEA.setLogStream(&logFile);
    testEA.setElitism(2);
    testEA.setDebugging(true);
    testEA.setPopulation(randomTrees);
    testEA.runGenerations(1000);

    cerr << "Final generation:\n";
    testEA.exportGenomes(cerr);

    logFile.close();

    return 0;
}
