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

PhyloTree decodePrefixNotation(string);

int main(int argc, const char *argv[])
{
    if (argc != 2) {
        cerr << "usage: " << argv[0] << " [fasta file]\n";
        return 1;
    }

    nodes = Fasta::readFastaFile(argv[1], true);

    string s;
    while (getline(cin, s)) {
        decodePrefixNotation(s);
        cout << "1.0\n";
    }

    return 0;
}

PhyloTree decodePrefixNotation(string s)
{
    std::istringstream inpStream(s);
    std::istringstream ss;

    string t;
    unsigned int n;

    // root node
    PhyloTreeNode* root = new PhyloTreeNode();
    stack<PhyloTreeNode*> nodeStack;
    nodeStack.push(root);

    // the first node should be an HTU (the root)
    assert(inpStream >> t && t == "h");

    while (inpStream >> t) {

        // HTU nodes are called “h#”
        if (t == "h") {
            cerr << "HTU\n";

            PhyloTreeNode* temp = new PhyloTreeNode();
            nodeStack.top()->addChild(temp, 1.0);
            nodeStack.push(temp);
        } else {
            ss.clear();
            ss.str(t);
            ss >> n;

            // if not starting with ‘h’ the token should always be a number
            assert(!ss.fail());

            cerr << "OTU #" << n << '\n';

            assert(n >= 0 && n < nodes.size());

            nodeStack.top()->addChild(nodes[n], 1.0);
        }

        while (nodeStack.size() > 0 && nodeStack.top()->numChildren() == 2) {
            cerr << "POP" << nodeStack.size() << "\n";
            nodeStack.pop();
        }
    }

    assert(nodeStack.size() == 0);

    PhyloTree tree(root, NULL);
    cout << tree.dot();

    //PhyloTree tree(root, NULL);
    //cout << tree.dot();

    return tree;
}
