#ifndef PHYLOTREE_H_INCLUDED
#define PHYLOTREE_H_INCLUDED

#include <string>
#include <vector>
#include <ostream>
#include <map>

#include "EvolutionModel.h"

using std::string;
using std::vector;
using std::ostream;
using std::map;

// super class for all phylo tree nodes
class PhyloTreeNode {
public:
    PhyloTreeNode();
    PhyloTreeNode(string name, string states);
    ~PhyloTreeNode();
    const string getName(bool printHTUs = true) const;
    unsigned int noStates() const { return m_nStates; }
    int height() const;
    void addChild(PhyloTreeNode*, double);
    void removeChild(PhyloTreeNode*);

    vector<vector<double> > likelihood(EvolutionModel*);

    static int count;

    string links() const;

    // returns string representation of tree in the newick format
    string newick() const;

    string name() const { return m_name; }

    // returns left/right child
    PhyloTreeNode* left() const { return m_left; }
    PhyloTreeNode* right() const { return m_right; }

    void setParent(PhyloTreeNode* parent) { m_parent = parent; }
    PhyloTreeNode* getParent() const { return m_parent; }

    // returns true if node has no parents
    bool isRoot() const { return m_parent == NULL; }

    // returns true if node has no children
    bool isLeaf() const { return m_left == NULL and m_right == NULL; }

    // returns number of children
    unsigned int numChildren() const;

    static string prefixRepresentation(PhyloTreeNode* node);

    const string getStates() const { return m_states; }

    void setNumStates(unsigned int n);

    friend ostream& operator<<(ostream& out, const PhyloTreeNode& n)
    {
        return (out << n.m_name);
    }

    void setLeftDist(double t) { m_leftDist = t; }
    void setRightDist(double t) { m_rightDist = t; }
    double getLeftDist() const { return m_leftDist; }
    double getRightDist() const { return m_rightDist; }

    PhyloTreeNode* findChild(string name);

    vector<PhyloTreeNode*> selfAndDescendants();

protected:

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

    // cache P and L values
    map<char, double> plCache;
    map<char, double> prCache;
    vector<vector<double> > llCache;
    vector<vector<double> > lrCache;

    // record the branch lengths the cache is for so the cache can be
    // invalidated on changes
    double cachedForL;
    double cachedForR;
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

    void setEvolutionModel(EvolutionModel* m) { m_evModel = m; }

    // return root node
    PhyloTreeNode* getRoot() const { return m_rootNode; }

    friend ostream& operator<<(ostream& out, const PhyloTree& t);

    static PhyloTree decodePrefixNotation(vector<PhyloTreeNode*> nodes, string s, EvolutionModel* evModel);

    void removeNode(string name);

private:
    PhyloTreeNode* m_rootNode;
    EvolutionModel *m_evModel;
};

#endif
