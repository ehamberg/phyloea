#include <fstream>
#include <cassert>
#include <iostream>
#include <set>
#include <sstream>

#include "Fasta.h"

using std::ifstream;
using std::stringstream;
using std::ios;
using std::cerr;
using std::set;

using namespace std;

vector<PhyloTreeNode*> Fasta::readFastaFile(string filename)
{
    vector<PhyloTreeNode*> nodes;

    ifstream fastaFile(filename.c_str());
    assert(fastaFile.is_open());

    vector<string> descriptions;
    vector<string> data;

    bool readingComment = false;
    int n = -1;
    string line;

    while (fastaFile >> line and !fastaFile.eof()) {
      if (line[0] == '>' or line[0] == ';') {
          if (!readingComment) {
              readingComment = true;
              n++;
              descriptions.push_back(string());
              data.push_back(string());
          }

          descriptions.back() += line.substr(1);
      } else {
          assert(n>=0); // we should have read a comment before the data
          readingComment = false;

          data.back() += line;
      }
    }
    fastaFile.close();

    assert(data.size() == descriptions.size());

    cerr << "Read " << data.size() << " sequences:\n";
    for (unsigned int i = 0; i < data.size(); i++) {
        stringstream ss;
        ss << i;
        cerr << "\t[" << i << "]" << ss.str() << " (length: " << data[i].length() << '\n';
        nodes.push_back(new PhyloTreeNode(ss.str(), data[i]));
    }

    return nodes;
}
