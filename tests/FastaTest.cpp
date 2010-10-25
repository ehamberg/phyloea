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
    vector<PhyloTreeNode*> nodes = Fasta::readFastaFile("tests/aligned.fasta");

    ASSERT_EQ(7, nodes.size());

    vector<PhyloTreeNode*>::iterator it;
    for (it = nodes.begin(); it != nodes.end(); it++) {
        delete *it;
    }
}
