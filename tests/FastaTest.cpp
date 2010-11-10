#include "gtest/gtest.h"
#include <vector>

#include "Fasta.h"
#include "PhyloTree.h"

class FastaTest : public ::testing::Test {
protected:
    FastaTest()
    {
    }

    virtual ~FastaTest()
    {
    }
};

// test a simple tree for which the likelihood is known
TEST_F(FastaTest, ReadFile) {
    vector<PhyloTreeNode*> nodes = Fasta::readFastaFile("tests/aligned.fasta", false);

    ASSERT_EQ(7, nodes.size());

    for (unsigned int i = 0; i < nodes.size(); i++) {
        ASSERT_EQ(18160, nodes[i]->noStates());
    }

    vector<PhyloTreeNode*>::iterator it;
    for (it = nodes.begin(); it != nodes.end(); it++) {
        delete *it;
    }
}

TEST_F(FastaTest, ReadInvalidFileDeathTest) {
    ASSERT_DEATH(Fasta::readFastaFile("tests/invalid.fasta"), "n>=0");
}
