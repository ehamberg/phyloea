#include "PhyloTree.h"
#include "EvolutionModel.h"
#include "Fasta.h"
#include "EASystem.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <stack>

#include "EAOperators.h"
#include "AMaxOperators.h"

using std::cout;
using std::cerr;
using std::cin;
using std::string;
using std::vector;
using std::stack;

vector<PhyloTreeNode*> nodes;
PhyloTree t;

int main(int argc, const char *argv[])
{
    if (argc != 2) {
        cerr << "usage: " << argv[0] << " [fasta file]\n";
        return 1;
    }

    nodes = Fasta::readFastaFile(argv[1]);

    string s;
    while (getline(cin, s)) {
        cout << PhyloTree::decodePrefixNotation(nodes, s, new Kimura(10.0)).logLikelihood() << '\n';
        //cout << PhyloTree::decodePrefixNotation(nodes, s, new Kimura(10.0)).dot() << '\n';
    }

    // clean-up
    for (unsigned int i = 0; i < nodes.size(); i++) {
        delete nodes[i];
    }
    nodes.clear();

    return 0;
}
