#include "TreeOperators.h"
#include "utils.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

vector<string> findSubTree(const vector<string>& t, unsigned int n)
{
    assert(n < t.size());
    vector<string> sub;

    int arity = 0;
    int i = n;

    sub.push_back(t.at(i));
    if (sub.front() == "h") {
        arity -= 2;
    }

    while (arity < 0) {
        // length
        sub.push_back(t.at(++i));

        // node
        sub.push_back(t.at(++i));

        if (sub.back() == "h") {
            arity -= 2;
        }
        ++arity;
    }

    return sub;
}

/** remove the given leaves from the tree.
 */
void pruneTree(vector<string>& tokens, const vector<string>& leaves)
{
    vector<string>::const_iterator it;
    for (it = leaves.begin(); it != leaves.end(); ++it) {
        for (unsigned int i = 0; i < tokens.size(); i+=2) {
            if (tokens.at(i) == "h") {
                continue;
            }

            if (tokens.at(i) == *it) { // we found the token
                // erase node
                //cerr << "erasing " << *it << endl;
                tokens.erase(tokens.begin()+i);

                // erase branch length
                //cerr << "erasing " << tokens.at(i-1) << endl;
                tokens.erase(tokens.begin()+(i-1));

                // find previous HTU node
                unsigned int j = i-2;
                while (tokens.at(j) != "h") --j;
                if (j > 0) {
                    // erase node
                    tokens.erase(tokens.begin()+j);

                    // record and erase branch length
                    tokens.erase(tokens.begin()+(j-1));

                } else if (j == 0) {
                    // there is only one node left on this side of the root
                    tokens.erase(tokens.begin(), tokens.begin()+2);
                }
                continue;
            }
        }
    }
}

void MutateTree::mutate(vector<string>& genomes)
{
    vector<string>::iterator it;
    for (it = genomes.begin(); it != genomes.end(); ++it) {
        if (randZeroToOne() > m_mutationRate) continue;

        // pick a random branch length
        unsigned int n = (rand()%(m_numNodes+1))*2+1;
        string s;

        //cerr << *it << endl;
        //cerr << "n = " << n << endl;

        unsigned int startPos = 0;
        unsigned int endPos;

        while (n-- > 0) {
            startPos = it->find_first_of('\t', startPos)+1;
        }

        endPos = it->find_first_of('\t', startPos+1);

        istringstream iss(it->substr(startPos, (endPos-startPos)));
        double len;
        iss >> len;

        // add [-0.1,0.1] to length
        len += -0.1+(randZeroToOne()*0.2);

        if (len < 0.000001) len = 0.000001;

        ostringstream oss;
        oss << len;

        cerr << "bef:\t" << *it << endl;
        it->replace(startPos, (endPos-startPos), oss.str());
        cerr << "mut:\t" << *it << endl;
        //cerr << *it << endl;
        //cerr << "==========" << endl;
    }
}

/** recombination operator for trees represented in preorder notation.
 *
 * does prune-delete-graft to recombine two trees
 */
vector<string> RecombineTree::produceOffspring(const string& p1, const string& p2)
{
    vector<string> children;
    if (randZeroToOne() >= m_recombProb) {
        cerr << "p1:\t" << p1 << endl;
        cerr << "p2:\t" << p2 << endl;
        vector<string> tokens1 = tokenize(p1);
        vector<string> tokens2 = tokenize(p2);

        // # of nodes in p1
        unsigned int n = ceil(tokens1.size()/2.0);

        // find a random subtree root from p1
        int subRoot = (rand()%(n-1)+1)*2;

        // find the p1 sub tree
        vector<string> subTree = findSubTree(tokens1, subRoot);

        // find the leaves of the p1 sub tree
        vector<string> subLeaves;
        vector<string>::const_iterator it;
        for (it = subTree.begin(); it < subTree.end(); it+=2) {
            if (*it != "h") {
                subLeaves.push_back(*it);
            }
        }

        // remove the leaves found in the p1 sub tree from p2
        pruneTree(tokens2, subLeaves);

        // find a random node on the tree
        int p = (rand()%int(ceil(tokens2.size()/2.0)))*2-1;
        if (p < 1) p = 0;

        if (p == 0) { // h x l y
            // prepend an HTU node
            subTree.insert(subTree.begin(), string("h"));
            subTree.insert(subTree.begin()+1, convertToString(randZeroToOne()));
            subTree.insert(subTree.end(), convertToString(randZeroToOne()));
        } else { // x h y l
            subTree.insert(subTree.begin(), convertToString(randZeroToOne()));
            subTree.insert(subTree.begin()+1, string("h"));
            subTree.insert(subTree.begin()+2, convertToString(randZeroToOne()));
        }

        // finally, add the sub tree to p2 at the point p
        tokens2.insert(tokens2.begin()+p, subTree.begin(), subTree.end());

        cerr << "off:\t" << detokenize(tokens2) << endl;
        children.push_back(detokenize(tokens2));
    } else {
        children.push_back(p1);
        children.push_back(p2);
    }

    return children;
}
