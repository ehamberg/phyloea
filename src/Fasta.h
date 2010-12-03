#ifndef FASTA_H_INCLUDED
#define FASTA_H_INCLUDED

#include <vector>
#include <string>
#include "PhyloTree.h"

using std::vector;
using std::string;

class Fasta {
public:
    // read a fasta file and return a vector of tree nodes containing the
    // nucleotide data from the file
    static vector<PhyloTreeNode*> readFastaFile(string filename, int maxLen = -1);
};

#endif
