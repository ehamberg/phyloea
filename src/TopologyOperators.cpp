#include "TopologyOperators.h"
#include "utils.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

vector<string> findSubTree2(const vector<string>& t, unsigned int n)
{
    assert(n <= t.size());
    vector<string> sub;

    int arity = 0;
    int i = n;

    sub.push_back(t.at(i));
    if (sub.front() == "h") {
        arity -= 2;
    }

    while (arity < 0) {
        // node
        sub.push_back(t.at(++i));

        if (sub.back() == "h") {
            arity -= 2;
        }
        ++arity;
    }

    return sub;
}

void pruneTree2(vector<string>& tokens, const vector<string>& leaves)
{
    for (auto it = leaves.cbegin(); it != leaves.cend(); ++it) {
        for (unsigned int i = 0; i < tokens.size(); i++) {
            if (tokens.at(i) == "h") {
                continue;
            }

            if (tokens.at(i) == *it) {
                // we found the token
                // erase node
                tokens.erase(tokens.begin()+i);

                // find previous HTU node
                unsigned int j = i-1;

                while (tokens.at(j) != "h") --j;
                if (j > 0) {
                    // erase node
                    tokens.erase(tokens.begin()+j);
                } else if (j == 0) {
                    // there is only one node left on this side of the root
                    tokens.erase(tokens.begin(), tokens.begin()+1);
                }
                break;
            }
        }
    }
}

// mutation: swapping two leaf nodes
void MutateTopology::mutate(vector<string>& genomes)
{
    unsigned int rand1, rand2;

    for (auto it = genomes.begin(); it != genomes.end(); ++it) {
        vector<string> tokens = tokenize(*it);
        rand1 = rand()%tokens.size();
        rand2 = rand()%tokens.size();

        // find the first leaf node (i.e. not an internal node)
        while (tokens.at(rand1) == "h") {
            ++rand1;
            rand1 = rand1%tokens.size();
        }

        while (tokens.at(rand2) == "h") {
            ++rand2;
            rand2 = rand2%tokens.size();
        }

        if (rand1 != rand2) {
            string temp = tokens.at(rand1);
            tokens[rand1] = tokens.at(rand2);
            tokens[rand2] = temp;
        }
        *it = detokenize(tokens);
    }
}

/** recombination operator for trees represented in preorder notation.
 *
 * does prune-delete-graft to recombine two trees
 */
vector<string> RecombineTopology::produceOffspring(const string& p1, const string& p2)
{
    vector<string> children;

    if (randZeroToOne() <= m_recombProb) {
        vector<string> tokens1 = tokenize(p1);
        vector<string> tokens2 = tokenize(p2);

        // # of nodes
        unsigned int n = tokens1.size();

        // find a random subtree root from p1
        int subRoot = 0;
        while (subRoot == 0) subRoot = rand()%n;

        //cerr << "subRoot: " << subRoot << "(" << tokens1.at(subRoot) << ")\n";

        // find the p1 sub tree
        vector<string> subTree = findSubTree2(tokens1, subRoot);

        // find the leaves of the p1 sub tree
        vector<string> subLeaves;
        for (auto it = subTree.cbegin(); it < subTree.cend(); ++it) {
            if (*it != "h") {
                subLeaves.push_back(*it);
            }
        }

        //// remove the leaves found in the p1 sub tree from p2
        pruneTree2(tokens2, subLeaves);

        // find a random insertion point
        int p = rand()%tokens2.size();

        // new internal node
        subTree.insert(subTree.begin(), string("h"));

        // finally, add the sub tree to p2 at the point p
        tokens2.insert(tokens2.begin()+p, subTree.begin(), subTree.end());

        children.push_back(detokenize(tokens2));

    } else {
        children.push_back(p1);
        children.push_back(p2);
    }

    return children;
}
