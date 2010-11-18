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
    PhyloTreeNode(PhyloTreeNode* parent);
    PhyloTreeNode(PhyloTreeNode* parent, string name, string states);
    ~PhyloTreeNode();
    const string getName(bool printAnons = true) const;
    unsigned int noStates() const { return m_nStates; }
    int height() const;
    void addChild(PhyloTreeNode*, double);

    vector<vector<double> > likelihood(EvolutionModel*);

    static int count;

    string links() const;

    // returns string representation of tree in the newick format
    string newick() const;

    string name() const { return m_name; }

    void setParent(PhyloTreeNode* parent) { m_parent = parent; }

    // returns true if node has no parents
    bool isRoot() const { return m_parent == NULL; }

    // returns number of children
    unsigned int numChildren() const;

    friend ostream& operator<<(ostream& out, const PhyloTreeNode& n)
    {
        return (out << n.m_name);
    }

private:
    string m_name; // species/taxon name
    string m_states; // observed states for each site

    // link to parent
    PhyloTreeNode* m_parent;

    // child nodes together with a distance
    PhyloTreeNode* m_left;
    PhyloTreeNode* m_right;
    double m_leftDist;
    double m_rightDist;

    vector<double> leafLikelihood(char n) const;

    vector<vector<double> > m_likelihoods;
    unsigned int m_nStates;
};


// convenience class for managing a tree of PhyloTreeNodes
class PhyloTree {
public:
    PhyloTree() { m_rootNode = 0; m_evModel = 0; }
    PhyloTree(PhyloTreeNode* r, EvolutionModel* m) : m_rootNode(r), m_evModel(m) {}
    ~PhyloTree();

    // returns the tree's height
    int height() const;

    // returns the tree's total likelihood
    double logLikelihood();

    // returns string representation of tree in graphviz dot format
    string dot() const;

    // returns string representation of tree in the newick format
    string newick() const;

    // build a random tree from the given nodes
    void buildRandomTree(vector<PhyloTreeNode*> leaves);

    // return root node
    PhyloTreeNode* getRoot() const { return m_rootNode; }

    friend ostream& operator<<(ostream& out, const PhyloTree& t);

private:
    PhyloTreeNode* m_rootNode;
    EvolutionModel *m_evModel;
};

#endif
