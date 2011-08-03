#include <fstream>
#include <cassert>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include "Fasta.h"

using std::ifstream;
using std::stringstream;
using std::ios;
using std::cerr;
using std::set;

using namespace std;

vector<PhyloTreeNode*> Fasta::readFastaFile(string filename, int maxLen)
{
    vector<PhyloTreeNode*> nodes;

    ifstream fastaFile(filename.c_str());
    assert(fastaFile.is_open());

    vector<string> descriptions;
    vector<string> data;

    bool readingComment = false;
    int n = -1;
    string line;

    while ( getline(fastaFile, line) and !fastaFile.eof()) {
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

          // remove space characters
          for (auto it = line.begin(); it != line.end(); ++it) {
              if (*it == ' ') line.erase(it);
          }
          data.back() += line;
      }
    }
    fastaFile.close();

    assert(data.size() == descriptions.size());

    //cerr << "FASTA: Read " << data.size() << " sequences:\n";
    for (unsigned int i = 0; i < data.size(); i++) {
        stringstream ss;
        ss << i;
        //cerr << "\t[" << i << "]" << ss.str() << " (length: " << data[i].length() << ")\n";

        // remove everything over maxlen
        if (maxLen != -1) {
            data[i] = data[i].substr(0,maxLen);
        }
        nodes.push_back(new PhyloTreeNode(ss.str(), data[i]));
    }

    return nodes;
}
