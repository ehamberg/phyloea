#include "PhyloTree.h"
#include "EvolutionModel.h"
#include "Fasta.h"
#include "EASystem.h"
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "EAOperators.h"
#include "AMaxOperators.h"

using std::cout;
using std::cerr;

int main(int argc, const char *argv[])
{
    srand(time(NULL));

    vector<PhyloTreeNode*> nodes = Fasta::readFastaFile("tests/aligned.fasta", true);
    PhyloTree t;
    t.buildRandomTree(nodes);
    cout << t.dot() << '\n';
    cerr << PhyloTreeNode::prefixRepresentation(t.getRoot()) << '\n';

    //vector<string> randomPop;
    //for (int i = 0; i < 10; i++) {
    //    string s;
    //    for (int j = 0; j < 20; j++) {
    //        s += '0';//(char)(rand()%2+(int)'0');
    //    }
    //    randomPop.push_back(s);
    //}

    //EASystem<string> testEA(new MutateString(0.01), new RecombineString(0.7), new RouletteWheelSelection<string>, new PipeFitnessFunc<string>);
    ////testEA.setLogStream(&cerr);
    //testEA.setElitism(4);
    //testEA.setDebugging(true);
    //testEA.setPopulation(randomPop);
    //testEA.runGenerations(1);

    //cerr << "Final generation:\n";
    //testEA.exportGenomes(cerr);

    return 0;
}
