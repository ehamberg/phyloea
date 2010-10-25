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
    PhyloTreeNode(string name, string states);
    ~PhyloTreeNode();
    const string getName() const { return name; }
    unsigned int noStates() const { return nStates; }
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
    string states; // observed states for each site

    // child nodes together with a distance
    PhyloTreeNode* left;
    PhyloTreeNode* right;
    double leftDist;
    double rightDist;

    vector<vector<double> > likelihoods;
    vector<double> leafLikelihood(char n) const;
    unsigned int nStates;
};


// convenience class for managing a tree of PhyloTreeNodes
class PhyloTree {
public:
    PhyloTree(PhyloTreeNode* r, EvolutionModel* m) : rootNode(r), evModel(m) {}
    ~PhyloTree();

    // returns the tree's height
    int height() const;

    // returns the tree's total likelihood
    double logLikelihood();

    // returns string representation of tree in graphviz dot format
    string dot() const;

    friend ostream& operator<<(ostream& out, const PhyloTree& t);

private:
    PhyloTreeNode* rootNode;
    EvolutionModel *evModel;
};

#endif
