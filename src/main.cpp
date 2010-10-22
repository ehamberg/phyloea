#include "PhyloTree.h"
#include "EvolutionModel.h"
#include <iostream>

using std::cout;

int main(int argc, const char *argv[])
{
    PhyloTreeNode *tree1 = new PhyloTreeNode();
    PhyloTreeNode *tree2 = new PhyloTreeNode();

    tree1->addChild(new PhyloTreeNode(), 1.0);
    tree1->addChild(new PhyloTreeNode(), 1.0);
    tree1->addChild(new PhyloTreeNode(), 1.0);

    tree2->addChild(tree1, 1.0);
    tree2->addChild(tree1, 1.0);
    tree2->addChild(tree1, 1.0);

    cout << tree1->height() << '\n';
    cout << tree2->height() << '\n';

    delete tree1;
    delete tree2;

    vector<char> states1;
    states1.push_back('A');
    PhyloTreeNode* n1 = new PhyloTreeNode("species 1", states1);

    vector<char> states2;
    states2.push_back('C');
    PhyloTreeNode* n2 = new PhyloTreeNode("species 2", states2);

    PhyloTreeNode y;
    y.addChild(n1, 0.5);
    y.addChild(n2, 0.5);
    Kimura *k = new Kimura(10);
    vector<vector<double> > foo = y.likelihood(k);
    cout << foo[0][0] << '\n';
    cout << foo[0][1] << '\n';
    cout << foo[0][2] << '\n';
    cout << foo[0][3] << '\n';

    cout << k->P('A', 'A', 0.5);

    return 0;
}
