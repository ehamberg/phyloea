#ifndef PHYLOTREE_H_INCLUDED
#define PHYLOTREE_H_INCLUDED

#include <string>
#include <vector>

#include "EvolutionModel.h"

using std::string;
using std::vector;
using std::pair;

// super class for all phylo tree nodes
class PhyloTreeNode {
public:
    PhyloTreeNode();
    PhyloTreeNode(string name, vector<char> states);
    const string getName() const { return name; }
    int noStates() const { return 1/*states.size()*/; }
    int height() const;
    void addChild(PhyloTreeNode*, double);

    vector<vector<double> > likelihood(EvolutionModel*);

    static int count;

private:
    string name; // species/taxon name
    vector<char> states; // observed states for each site
    vector<pair<PhyloTreeNode*, double> > children; // child nodes together with a distance
    vector<vector<double> > likelihoods;
    vector<double> leafLikelihood(char n) const;
};


// convenience class for managing a tree of PhyloTreeNodes
class PhyloTree {
    int height() const;

private:
    PhyloTreeNode* rootNode;
};

#endif
