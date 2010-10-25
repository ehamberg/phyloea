#include "gtest/gtest.h"

#include <vector>
#include <cmath>

#include "PhyloTree.h"
#include "EvolutionModel.h"

class PhyloTreeTest : public ::testing::Test {
protected:
    PhyloTreeTest() {
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

        PhyloTreeNode* y = new PhyloTreeNode();
        y->addChild(n1, 0.5);
        y->addChild(n2, 0.5);

        PhyloTreeNode* w = new PhyloTreeNode();
        w->addChild(n4, 0.5);
        w->addChild(n5, 0.5);

        PhyloTreeNode* z = new PhyloTreeNode();
        z->addChild(n3, 0.5);
        z->addChild(w, 0.5);

        PhyloTreeNode* x = new PhyloTreeNode();
        x->addChild(y, 0.5);
        x->addChild(z, 0.5);

        Kimura *k = new Kimura(10);
        t = new PhyloTree(x, k);

        tree1 = new PhyloTreeNode();
        tree2 = new PhyloTreeNode();

        tree1->addChild(new PhyloTreeNode(), 1.0);
        tree1->addChild(new PhyloTreeNode(), 1.0);
        tree1->addChild(new PhyloTreeNode(), 1.0);

        tree2->addChild(tree1, 1.0);
        tree2->addChild(tree1, 1.0);
        tree2->addChild(tree1, 1.0);
    }

    virtual ~PhyloTreeTest() {
        delete t;
        delete tree1;
        delete tree2;
    }

    PhyloTree* t;
    PhyloTreeNode *tree1;
    PhyloTreeNode *tree2;


};

// test a simple tree for which the likelihood is known
TEST_F(PhyloTreeTest, SimpleTreeLikelihood) {
    ASSERT_FLOAT_EQ(0.0000294480138762, pow(10.0, t->logLikelihood()));
}

TEST_F(PhyloTreeTest, TreeHeight) {
    ASSERT_EQ(2, tree1->height());
    ASSERT_EQ(3, tree2->height());
}
