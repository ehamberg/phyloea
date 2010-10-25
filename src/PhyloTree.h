#ifndef PHYLOTREE_H_INCLUDED
#define PHYLOTREE_H_INCLUDED

#include <string>
#include <vector>
#include <ostream>

#include "EvolutionModel.h"

using std::string;
using std::vector;
using std::pair;
using std::ostream;

// super class for all phylo tree nodes
class PhyloTreeNode {
public:
    PhyloTreeNode();
    PhyloTreeNode(string name, vector<char> states);
    ~PhyloTreeNode();
    const string getName() const { return name; }
    int noStates() const { return 1/*states.size()*/; }
    int height() const;
    void addChild(PhyloTreeNode*, double);

    vector<vector<double> > likelihood(EvolutionModel*);

    static int count;

    string links() const;

    friend ostream& operator<<(ostream& out, const PhyloTreeNode& n)
    {
        return (out << n.name);
    }

private:
    string name; // species/taxon name
    vector<char> states; // observed states for each site
    vector<pair<PhyloTreeNode*, double> > children; // child nodes together with a distance
    vector<vector<double> > likelihoods;
    vector<double> leafLikelihood(char n) const;
};


// convenience class for managing a tree of PhyloTreeNodes
class PhyloTree {
public:
    PhyloTree(PhyloTreeNode* r, EvolutionModel* m) : rootNode(r), evModel(m) {}
    ~PhyloTree();

    // returns the tree's height
    int height() const;

    // returns the tree's total likelihood
    double likelihood();

    // returns string representation of tree in graphviz dot format
    string dot() const;

    friend ostream& operator<<(ostream& out, const PhyloTree& t);

private:
    PhyloTreeNode* rootNode;
    EvolutionModel *evModel;
};

#endif
