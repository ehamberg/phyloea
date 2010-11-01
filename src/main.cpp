#include "PhyloTree.h"
#include "EvolutionModel.h"
#include "Fasta.h"
#include <iostream>
#include <cmath>

using std::cout;

int main(int argc, const char *argv[])
{
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

    return 0;
}
