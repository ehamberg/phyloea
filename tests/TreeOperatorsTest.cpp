#include "gtest/gtest.h"

#include <vector>
#include <cmath>
#include <iostream>

#include "PhyloTree.h"
#include "EvolutionModel.h"
#include "TreeOperators.h"
#include "Fasta.h"
#include "utils.h"

class TreeOperatorsTest : public ::testing::Test {
protected:
    TreeOperatorsTest() {
    }

    virtual ~TreeOperatorsTest() {
    }
};

// make sure branch lenghts are kept for nodes moved “up” in a tree as it is
// pruned
TEST_F(TreeOperatorsTest, PruneLengths1) {
    vector<string> tokens;
    vector<string> leaves;

    tokens.resize(13);
    string toks[] = {"h", "7", "A", "8", "h", "9", "h", "11", "B", "12", "C",
        "10", "D"};
    copy (toks, toks+13, tokens.begin());

    leaves.push_back("C");
    leaves.push_back("D");

    pruneTree(tokens, leaves);

    ASSERT_EQ("h\t7\tA\t11\tB", detokenize(tokens));
}

// make sure branch lenghts are kept for nodes moved “up” in a tree as it is
// pruned
TEST_F(TreeOperatorsTest, PruneLengths2) {
    vector<string> tokens;
    vector<string> leaves;

    tokens.resize(13);
    string toks[] = {"h", "7", "A", "8", "h", "9", "h", "11", "B", "12", "C",
        "10", "D"};
    copy (toks, toks+13, tokens.begin());

    leaves.push_back("A");
    leaves.push_back("D");

    pruneTree(tokens, leaves);

    ASSERT_EQ("h\t11\tB\t12\tC", detokenize(tokens));
}

TEST_F(TreeOperatorsTest, PruneLengths3) {
    vector<string> tokens;
    vector<string> leaves;

    tokens.resize(5);
    string toks[] = {"h", "7", "A", "8", "B"};
    copy (toks, toks+5, tokens.begin());

    leaves.push_back("A");

    pruneTree(tokens, leaves);

    ASSERT_EQ("B", detokenize(tokens));
}
