#include <fstream>
#include <cassert>
#include <iostream>
#include <set>

#include "Fasta.h"

using std::ifstream;
using std::ios;
using std::cerr;
using std::set;

using namespace std;

vector<PhyloTreeNode*> Fasta::readFastaFile(string filename, bool removeGaps)
{
    vector<PhyloTreeNode*> nodes;

    ifstream fastaFile(filename.c_str());
    assert(fastaFile.is_open());

    vector<string> descriptions;
    vector<string> data;

    bool readingComment = false;
    int n = -1;
    string line;

    while (fastaFile >> line && !fastaFile.eof()) {
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

    if (removeGaps) {
        set<unsigned int> gapSites;

        // record all gap sites ...
        vector<string>::iterator it;
        for (it = data.begin(); it != data.end(); ++it) {
            for (unsigned int i = 0; i < it->length(); i++) {
                if (it->at(i) == '-') {
                    gapSites.insert(i);
                }
            }
        }

        cerr << "Removing " << gapSites.size() << " gaps...\n";

        // ... then remove those sites from all sequences
        set<unsigned int>::iterator jt;
        for (it = data.begin(); it != data.end(); ++it) {
            unsigned int i = 0;
            for (jt = gapSites.begin(); jt != gapSites.end(); ++jt) {
                cerr << (*jt+i) << '\n';
                it->erase(*jt+(i--), 1);
            }
        }
    }

    cerr << "Read " << data.size() << " sequences:\n";
    for (unsigned int i = 0; i < data.size(); i++) {
        cerr << '\t' << descriptions[i] << " (length: " << data[i].length() << '\n';
    }


    for (unsigned int i = 0; i < data.size(); i++) {
        nodes.push_back(new PhyloTreeNode(descriptions[i], data[i]));
    }

    return nodes;
}
