#include "gtest/gtest.h"

#include <vector>
#include <iostream>

#include "../src/PhyloTree.h"
#include "../src/EvolutionModel.h"

namespace {

class PhyloTreeTest : public ::testing::Test {
protected:
    PhyloTreeTest() {
        std::vector<char> states1;
        states1.push_back('A');
        PhyloTreeNode* n1 = new PhyloTreeNode("species 1", states1);

        std::vector<char> states2;
        states2.push_back('C');
        PhyloTreeNode* n2 = new PhyloTreeNode("species 2", states2);

        y = new PhyloTreeNode();
        y->addChild(n1, 0.5);
        y->addChild(n2, 0.5);
    }

    virtual ~PhyloTreeTest() {
        delete y;
    }

    PhyloTreeNode* y;
};

// test a simple tree for which the likelihood is known
TEST_F(PhyloTreeTest, SimpleTreeLikelihood) {
    std::vector<double> likelihoods = y->likelihood(new Kimura(10))[0];
    double sum = 0;
    for (int i = 0; i < 4; i++) {
        sum += likelihoods[i];
    }
    ASSERT_FLOAT_EQ(0.041561770481204886, sum);
}

}  // namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
