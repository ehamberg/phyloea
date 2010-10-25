#ifndef FASTA_H_INCLUDED
#define FASTA_H_INCLUDED

#include <vector>
#include <string>
#include "PhyloTree.h"

using std::vector;
using std::string;

class Fasta {
public:
    static vector<PhyloTreeNode*> readFastaFile(string filename);
};

#endif
