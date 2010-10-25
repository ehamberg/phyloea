#include "PhyloTree.h"
#include "EvolutionModel.h"
#include <iostream>

using std::cout;

int main(int argc, const char *argv[])
{
    vector<char> states;
    states.push_back('A');
    PhyloTreeNode* n1 = new PhyloTreeNode("species 1", states);

    states[0] = 'C';
    PhyloTreeNode* n2 = new PhyloTreeNode("species 2", states);

    states[0] = 'C';
    PhyloTreeNode* n3 = new PhyloTreeNode("species 3", states);

    states[0] = 'C';
    PhyloTreeNode* n4 = new PhyloTreeNode("species 4", states);

    states[0] = 'G';
    PhyloTreeNode* n5 = new PhyloTreeNode("species 5", states);

    PhyloTreeNode y;
    y.addChild(n1, 0.5);
    y.addChild(n2, 0.5);

    PhyloTreeNode w;
    w.addChild(n4, 0.5);
    w.addChild(n5, 0.5);

    PhyloTreeNode z;
    z.addChild(n3, 0.5);
    z.addChild(&w, 0.5);

    PhyloTreeNode x;
    x.addChild(&y, 0.5);
    x.addChild(&z, 0.5);

    Kimura *k = new Kimura(10);
    PhyloTree t(&x, k);

    cout << t.likelihood() << '\n';
    cout << t.dot() << '\n';

    delete k;

    return 0;
}
