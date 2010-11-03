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

int main(int argc, const char *argv[])
{
    srand(time(NULL));

    PhyloTreeNode* n1 = new PhyloTreeNode("species 1", "AA");

    PhyloTreeNode* n2 = new PhyloTreeNode("species 2", "CC");

    PhyloTreeNode* n3 = new PhyloTreeNode("species 3", "CC");

    PhyloTreeNode* n4 = new PhyloTreeNode("species 4", "CC");

    PhyloTreeNode* n5 = new PhyloTreeNode("species 5", "GG");

    PhyloTreeNode* y = new PhyloTreeNode;
    y->addChild(n1, 0.5);
    y->addChild(n2, 0.5);

    PhyloTreeNode* w = new PhyloTreeNode;
    w->addChild(n4, 0.5);
    w->addChild(n5, 0.5);

    PhyloTreeNode* z = new PhyloTreeNode;
    z->addChild(n3, 0.5);
    z->addChild(w, 0.5);

    PhyloTreeNode* x = new PhyloTreeNode;
    x->addChild(y, 0.5);
    x->addChild(z, 0.5);

    PhyloTree t(x, new Kimura(10));

    cout << pow(10, t.logLikelihood()) << '\n';
    cout << t.dot() << '\n';

    vector<PhyloTreeNode*> nodes = Fasta::readFastaFile("tests/aligned.fasta");

    vector<PhyloTreeNode*>::iterator it;
    for (it = nodes.begin(); it != nodes.end(); ++it) {
        delete *it;
    }

    vector<string> randomPop;
    for (int i = 0; i < 4; i++) {
        string s;
        for (int j = 0; j < 10; j++) {
            s += (char)(rand()%5+(int)'A');
        }
        randomPop.push_back(s);
    }

    cout << "ea\n";
    EASystem<string> testEA(new MutateString(0.1), new RecombineString(0.7), new RouletteWheelSelection<string>);
    testEA.setPopulation(randomPop);
    testEA.runGenerations(10);
    testEA.exportGenomes(cout);

    return 0;
}
