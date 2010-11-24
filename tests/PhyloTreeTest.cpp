#include "gtest/gtest.h"

#include <vector>
#include <cmath>
#include <iostream>

#include "PhyloTree.h"
#include "EvolutionModel.h"
#include "Fasta.h"

class PhyloTreeTest : public ::testing::Test {
protected:
    PhyloTreeTest() {
        PhyloTreeNode* n1 = new PhyloTreeNode("s1", "A");

        PhyloTreeNode* n2 = new PhyloTreeNode("s2", "C");

        PhyloTreeNode* n3 = new PhyloTreeNode("s3", "C");

        PhyloTreeNode* n4 = new PhyloTreeNode("s4", "C");

        PhyloTreeNode* n5 = new PhyloTreeNode("s5", "G");

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
    ASSERT_FLOAT_EQ(0.00000011503131, pow(10.0, t->logLikelihood()));
}

TEST_F(PhyloTreeTest, NewickExportTest) {
    ASSERT_EQ("((s1:0.5,s2:0.5):0.5,(s3:0.5,(s4:0.5,s5:0.5):0.5):0.5);", t->newick());
}

TEST_F(PhyloTreeTest, TreeHeight) {
    ASSERT_EQ(2, tree1->height());
    ASSERT_EQ(3, tree2->height());
}

TEST_F(PhyloTreeTest, NumStates) {
    ASSERT_EQ(1, t->getRoot()->noStates());
}

TEST_F(PhyloTreeTest, RandomTree) {
    vector<PhyloTreeNode*> leaves;
    unsigned int n = 30;

    for (unsigned int i = 0; i < n; i++) {
        leaves.push_back(new PhyloTreeNode());
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

TEST_F(PhyloTreeTest, PulleyPrinciple) {
    double t1 = rand()/(double)RAND_MAX;
    double t2 = rand()/(double)RAND_MAX;
    PhyloTreeNode* r1 = new PhyloTreeNode();
    PhyloTreeNode* h1 = new PhyloTreeNode();

    h1->addChild(new PhyloTreeNode("s1", "ACGT"), rand()/(double)RAND_MAX);
    h1->addChild(new PhyloTreeNode("s2", "GCCA"), rand()/(double)RAND_MAX);

    r1->addChild(new PhyloTreeNode("s3", "GGGA"), t1);
    r1->addChild(h1, t2);

    PhyloTree tr1(r1, new Kimura(10.0));

    double lh1 = tr1.logLikelihood();
    ASSERT_LT(tr1.logLikelihood(), 0.0);

    for (unsigned int i = 0; i < 1000; i++) {
        double t_rand = rand()/(double)RAND_MAX;
        r1->setLeftDist(t_rand);
        r1->setLeftDist((t1+t2)-t_rand);
        ASSERT_FLOAT_EQ(tr1.logLikelihood(), lh1);
    }
}

TEST_F(PhyloTreeTest, PrefixCoding) {
    vector<PhyloTreeNode*> nodes1 = Fasta::readFastaFile("tests/aligned.fasta");
    PhyloTree t;
    t.setEvolutionModel(new Kimura(10.0));
    t.buildRandomTree(nodes1);

    double lh1 = t.logLikelihood();
    ASSERT_LT(lh1, 0.0);

    string prefixCoded = PhyloTreeNode::prefixRepresentation(t.getRoot());

    vector<PhyloTreeNode*> nodes2 = Fasta::readFastaFile("tests/aligned.fasta");
    PhyloTree t2 = PhyloTree::decodePrefixNotation(nodes2, prefixCoded, new Kimura(10.0));

    double lh2 = t.logLikelihood();

    ASSERT_DOUBLE_EQ(lh1, lh2);
}

