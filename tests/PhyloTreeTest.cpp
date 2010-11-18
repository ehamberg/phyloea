#include "gtest/gtest.h"

#include <vector>
#include <cmath>

#include "PhyloTree.h"
#include "EvolutionModel.h"

class PhyloTreeTest : public ::testing::Test {
protected:
    PhyloTreeTest() {
        PhyloTreeNode* n1 = new PhyloTreeNode(NULL, "s1", "A");

        PhyloTreeNode* n2 = new PhyloTreeNode(NULL, "s2", "C");

        PhyloTreeNode* n3 = new PhyloTreeNode(NULL, "s3", "C");

        PhyloTreeNode* n4 = new PhyloTreeNode(NULL, "s4", "C");

        PhyloTreeNode* n5 = new PhyloTreeNode(NULL, "s5", "G");

        PhyloTreeNode* y = new PhyloTreeNode(NULL);
        y->addChild(n1, 0.5);
        y->addChild(n2, 0.5);

        PhyloTreeNode* w = new PhyloTreeNode(NULL);
        w->addChild(n4, 0.5);
        w->addChild(n5, 0.5);

        PhyloTreeNode* z = new PhyloTreeNode(NULL);
        z->addChild(n3, 0.5);
        z->addChild(w, 0.5);

        PhyloTreeNode* x = new PhyloTreeNode(NULL);
        x->addChild(y, 0.5);
        x->addChild(z, 0.5);

        Kimura *k = new Kimura(10);
        t = new PhyloTree(x, k);

        tree1 = new PhyloTreeNode(NULL);
        tree2 = new PhyloTreeNode(NULL);

        tree1->addChild(new PhyloTreeNode(tree1), 1.0);
        tree1->addChild(new PhyloTreeNode(tree1), 1.0);

        tree2->addChild(tree1, 1.0);
    }

    virtual ~PhyloTreeTest() {
        delete t;
        delete tree2;
    }

    PhyloTree* t;
    PhyloTreeNode* tree1;
    PhyloTreeNode* tree2;
};

// test a simple tree for which the likelihood is known
TEST_F(PhyloTreeTest, SimpleTreeLikelihood) {
    ASSERT_FLOAT_EQ(0.0000294480138762, pow(10.0, t->logLikelihood()));
}

TEST_F(PhyloTreeTest, NewickExportTest) {
    ASSERT_EQ("((s1:0.5,s2:0.5):0.5,(s3:0.5,(s4:0.5,s5:0.5):0.5):0.5);", t->newick());
}

TEST_F(PhyloTreeTest, TreeHeight) {
    ASSERT_EQ(2, tree1->height());
    ASSERT_EQ(3, tree2->height());
}

TEST_F(PhyloTreeTest, RandomTree) {
    vector<PhyloTreeNode*> leaves;
    unsigned int n = 30;

    for (unsigned int i = 0; i < n; i++) {
        leaves.push_back(new PhyloTreeNode(NULL));
    }

    PhyloTree t;
    t.buildRandomTree(leaves);

    ASSERT_LT(log2(n), t.height()); // |tree| must be ≥ log₂(n)
    ASSERT_GE(n+1, t.height());     // |tree| must be ≤ n+1
    ASSERT_TRUE(t.getRoot()->isRoot()); // there should be a parentless root node

    // all leaves should now have a parent
    for (unsigned int i = 0; i < n; i++) {
        ASSERT_TRUE(!leaves.at(i)->isRoot());
    }
}
