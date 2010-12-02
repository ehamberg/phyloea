#include "TreeOperators.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

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

        len += -0.1+(rand()/(double)RAND_MAX*0.2);

        if (len < 0.0) len = 0.0;

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
        std::cerr << "lmao\t" << p1 << ", " << p2 << std::endl;
        assert(false);
    } else {
      children.push_back(p1);
      children.push_back(p2);
    }

    return children;
}
