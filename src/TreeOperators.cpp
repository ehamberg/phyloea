#include "TreeOperators.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <stdexcept>

using namespace std;

// FIXME: move to utils.cpp
vector<string> tokenize(const string& s)
{
    istringstream iss(s);
    vector<string> tokens;
    string temp;

    while (iss >> temp) tokens.push_back(temp);

    return tokens;
}

class BadConversion : public runtime_error {
public:
    BadConversion(string const& s) : runtime_error(s) {}
};

inline double convertToDouble(string const& s)
{
    std::istringstream i(s);
    double x;
    if (!(i >> x))
        throw BadConversion("convertToDouble(\"" + s + "\")");
    return x;
}

string detokenize(const vector<string>& tokens, string sep = string("\t"))
{
    string s;

    vector<string>::const_iterator it;
    for (it = tokens.begin(); it != tokens.end(); ++it) {
        s += *it;
        if (it != tokens.end()-1) {
            s += sep;
        }
    }

    return s;
}

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

void pruneTree(vector<string>& tokens, const vector<string>& leaves)
{
    vector<string>::const_iterator it;
    for (it = leaves.begin(); it != leaves.end(); ++it) {
        for (unsigned int i = 0; i < tokens.size(); i+=2) {
            if (tokens.at(i) == "h") {
                continue;
            }

            if (tokens.at(i) == *it) {
                // we found the token
                // erase node
                tokens.erase(tokens.begin()+i);
                // erase branch length
                tokens.erase(tokens.begin()+(i-1));

                // find previous HTU node
                unsigned int j = i-2;
                while (tokens.at(j) != "h") --j;
                if (j > 0) {
                    // erase node
                    tokens.erase(tokens.begin()+j);

                    // record and erase branch length
                    string len = tokens.at(j-1);
                    tokens.erase(tokens.begin()+(j-1));

                    double newLen = convertToDouble(len)+convertToDouble(tokens.at(j-1));
                    std::ostringstream s;
                    s << newLen;
                    tokens[j-1] = s.str();
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
        len += -0.1+(rand()/(double)RAND_MAX*0.2);

        if (len < 0.000001) len = 0.000001;

        ostringstream oss;
        oss << len;

        it->replace(startPos, (endPos-startPos), oss.str());
        //cerr << *it << endl;
        //cerr << "==========" << endl;
    }
}

vector<string> RecombineTree::produceOffspring(const string& p1, const string& p2)
{
    vector<string> children;
    if ((double)rand()/(double)RAND_MAX >= m_recombProb) {
        vector<string> tokens1 = tokenize(p1);
        vector<string> tokens2 = tokenize(p2);

        // # of nodes
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

        int p = (rand()%int(ceil(tokens2.size()/2.0)))*2-1;
        if (p < 1) p = 0;

        if (p == 0) { // h x l y
            // prepend an HTU node
            subTree.insert(subTree.begin(), string("h"));
            subTree.insert(subTree.begin()+1, string("0.666"));
            subTree.insert(subTree.end(), string("0.666"));
        } else { // x h y l
            subTree.insert(subTree.begin(), string("0.666"));
            subTree.insert(subTree.begin()+1, string("h"));
            subTree.insert(subTree.begin()+2, string("0.666"));
        }
        // finally, add the sub tree to p2 at the point p
        tokens2.insert(tokens2.begin()+p, subTree.begin(), subTree.end());

        children.push_back(detokenize(tokens2));
    } else {
        children.push_back(p1);
        children.push_back(p2);
    }

    return children;
}
