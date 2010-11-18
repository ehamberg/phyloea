#include "PhyloTree.h"
#include "EvolutionModel.h"
#include "Fasta.h"
#include "EASystem.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "EAOperators.h"
#include "AMaxOperators.h"

using std::cout;
using std::cerr;
using std::cin;
using std::string;
using std::vector;

vector<PhyloTreeNode*> nodes;
PhyloTree t;

PhyloTree decodePrefixNotation(string);

int main(int argc, const char *argv[])
{
    if (argc != 2) {
        cerr << "usage: " << argv[0] << " [fasta file]\n";
        return 1;
    }

    nodes = Fasta::readFastaFile(argv[1], true);
    t.buildRandomTree(nodes);

    string s;
    while (getline(cin, s)) {
        decodePrefixNotation(s);
        cout << "1.0\n";
    }

    return 0;
}

vector<string> tokenize_str(const string & str)
{
    char delims='\t';
    // Skip delims at beginning, find start of first token
    string::size_type lastPos = str.find_first_not_of(delims, 0);
    // Find next delimiter @ end of token
    string::size_type pos     = str.find_first_of(delims, lastPos);

    // output vector
    vector<string> tokens;

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delims.  Note the "not_of". this is beginning of token
        lastPos = str.find_first_not_of(delims, pos);
        // Find next delimiter at end of token.
        pos     = str.find_first_of(delims, lastPos);
    }

    return tokens;
}

PhyloTree decodePrefixNotation(string s)
{
    //vector<string> tokens = tokenize_str(s);

    std::istringstream inpStream(s);

    string t;
    int i;
    while (inpStream >> t) {
        cout << t << '\n';
    }

    return PhyloTree();
}
