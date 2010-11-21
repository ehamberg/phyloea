#include "PhyloTree.h"

#include <sstream>
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>

using std::vector;
using std::string;
using std::ostringstream;
using std::cerr;
using std::map;

#define MAX(x,y) (x>y?x:y)

int PhyloTreeNode::count = 0;

PhyloTreeNode::PhyloTreeNode()
{
    m_parent = NULL;

    // not given a name, so this is an HTU. generate a name for this node of
    // the form “h#”
    ostringstream out;
    out << PhyloTreeNode::count++;
    out.str();

    m_name = 'n'+out.str();
    m_nStates = 0;

    m_left = m_right = NULL;
}

PhyloTreeNode::PhyloTreeNode(string name, string states)
{
    m_parent = NULL;
    m_name = name;
    m_states = states;
    m_nStates = m_states.size();

    m_left = m_right = NULL;
}

PhyloTreeNode::~PhyloTreeNode()
{
    // delete children
    delete m_left;
    delete m_right;
}

void PhyloTreeNode::addChild(PhyloTreeNode* child, double distance) {
    assert(!m_left or !m_right); // only two child nodes are allowed

    child->setParent(this);

    if (m_nStates == 0) {
        m_nStates = child->noStates();
    } else {
        assert(m_nStates == child->noStates());
    }

    if (!m_left) {
        m_left = child;
        m_leftDist = distance;
    } else if (!m_right) {
        m_right = child;
        m_rightDist = distance;
    }
}

const string PhyloTreeNode::getName(bool printHTUs) const
{
    if (!printHTUs && numChildren() > 0) {
        return string();
    } else {
        return m_name;
    }
}

int PhyloTreeNode::height() const
{
    if (!m_left and !m_right) {
        return 1;
    }

    if (!m_right) {
        return 1+m_left->height();
    }

    return 1+MAX(m_left->height(), m_right->height());
}

string PhyloTreeNode::prefixRepresentation(PhyloTreeNode* node)
{
    if (node->numChildren() == 0) {
        return node->getName();
    }

    assert(node->numChildren() == 2);
    return node->getName()
        + '\t' + prefixRepresentation(node->left())
        + '\t' + prefixRepresentation(node->right());
}


int PhyloTree::height() const
{
    return m_rootNode->height();
}

vector<double> PhyloTreeNode::leafLikelihood(char n) const {
    vector<double> m_likelihoods;
    for (int i = 0; i < 4; i++) {
        m_likelihoods.push_back(0.0);
    }

    switch (n) {
    case 'A':
        m_likelihoods[0] = 1.0;
        break;
    case 'C':
        m_likelihoods[1] = 1.0;
        break;
    case 'G':
        m_likelihoods[2] = 1.0;
        break;
    case 'T':
        m_likelihoods[3] = 1.0;
        break;
    }

    return m_likelihoods;
}

vector<vector<double> > PhyloTreeNode::likelihood(EvolutionModel* eModel)
{
    //  base case: we have already found likelihoods for the node
    if (!m_likelihoods.empty()) {
        return m_likelihoods;
    }

    if (!m_left and !m_right) { // leaf node
        for (unsigned int i = 0; i < noStates(); i++) {
            m_likelihoods.push_back(leafLikelihood(m_states[i]));
        }

        return m_likelihoods;
    }

    // check if the cached values are for a different branch length, if so,
    // invalidate the caches
    if (cachedForL != m_leftDist) {
        plCache.clear();
        llCache.clear();
    }
    if (cachedForR != m_rightDist) {
        prCache.clear();
        lrCache.clear();
    }

    map<char, double>::iterator it;

    double xs, ys;
    char nucleotides[4] = {'A', 'C', 'G', 'T' };

    cerr << "computing for " << noStates() << " states...\n";
    for (unsigned int i = 0; i < noStates(); i++) {
        vector<double> l;

        if (i%1000==0) {
            cerr << '\t' << i << '\n';
        }
        for (int s = 0; s < 4; s++) {
            xs = ys = 0.0;
            for (int x = 0; x < 4; x++) {
                // p1
                double p1;
                it = plCache.find(nucleotides[s]^nucleotides[x]);
                if (it != plCache.end()) {
                    p1 = plCache[nucleotides[s]^nucleotides[x]];
                } else {
                    // not in cache
                    p1 = eModel->P(nucleotides[s], nucleotides[x], m_leftDist);
                    plCache[nucleotides[s]^nucleotides[x]] = p1;
                }

                // p2
                double p2;
                it = prCache.find(nucleotides[s]^nucleotides[x]);
                if (it != prCache.end()) {
                    p2 = prCache[nucleotides[s]^nucleotides[x]];
                } else {
                    // not in cache
                    p2 = eModel->P(nucleotides[s], nucleotides[x], m_rightDist);
                    prCache[nucleotides[s]^nucleotides[x]] = p2;
                }

                if (llCache.empty()) {
                    assert(lrCache.empty());
                    cachedForL = m_leftDist;
                    cachedForR = m_rightDist;
                    llCache = m_left->likelihood(eModel);
                    lrCache = m_right->likelihood(eModel);
                }
                double l1 = llCache[i][x];
                double l2 = lrCache[i][x];
                xs += p1*l1;
                ys += p2*l2;
            }
            l.push_back(eModel->prior(nucleotides[s])*xs*ys);
        }

        m_likelihoods.push_back(l);
    }
    cerr << "... done\n";

    return m_likelihoods;
}

// print all links from this node to child nodes in graphviz format
string PhyloTreeNode::links() const
{
    ostringstream out;

    if (m_left) {
        out << "\t\"" << getName() << "\" -> \"" << m_left->getName()
            << "\" [label=\"" << m_leftDist << "\"];\n";
        out << m_left->links();
    }

    if (m_right) {
        out << "\t\"" << getName() << "\" -> \"" << m_right->getName()
            << "\" [label=\"" << m_leftDist << "\"];\n";
        out << m_right->links();
    }

    return out.str();
}

string PhyloTreeNode::newick() const
{
    char s[256];

    if (!m_left or !m_right) {
        assert(!m_left and !m_right);
        return getName(false);
    } else {
        std::sprintf(s, "(%s:%.1f,%s:%.1f)%s", m_left->newick().c_str(),
                m_leftDist, m_right->newick().c_str(), m_rightDist,
                getName(false).c_str());
    }

    return string(s);
}
unsigned int PhyloTreeNode::numChildren() const
{
    unsigned int c = 0;
    if (m_left != NULL) c++;
    if (m_right != NULL) c++;

    return c;
}

ostream& operator<<(ostream& out, const PhyloTree& t)
{
    out << t.m_rootNode->links();

    return out;
}

string PhyloTree::dot() const
{
    return "digraph {\n" + m_rootNode->links() + "}\n";
}

string PhyloTree::newick() const
{
    PhyloTreeNode* t = m_rootNode;

    string newick = t->newick()/*+t->name()*/+';';
    std::replace(newick.begin(), newick.end(), ' ', '_');
    //"(((s4:05,s5:0.5)node_1:0.5,s3:0.5)node_2:0.5,(s1:0.5,s2:0.5)node_0:0.5)node_3;"

    return newick;
}

double PhyloTree::logLikelihood()
{
    double ret = 0.0;

    vector<vector<double> > siteLikelihoods = m_rootNode->likelihood(m_evModel);

    vector<vector<double> >::const_iterator it;
    vector<double>::const_iterator jt;
    for (it = siteLikelihoods.begin(); it != siteLikelihoods.end(); ++it) {

        double siteSum = 0.0;
        for (jt = (*it).begin(); jt != (*it).end(); ++jt) {
            siteSum += (*jt);
        }

        ret += log10(siteSum/4.0);
    }

    return ret;
}

PhyloTree::~PhyloTree()
{
    delete m_evModel;
    delete m_rootNode;
}

void PhyloTree::buildRandomTree(vector<PhyloTreeNode*> leaves)
{
    unsigned int len = leaves.size();

    // shuffle vector of leaves
    random_shuffle(leaves.begin(), leaves.end());

    // make a set of n-2 internal nodes + 1 root node ...
    vector<PhyloTreeNode*> internal;
    for (unsigned int i = 0; i < len-1; i++) {
        internal.push_back(new PhyloTreeNode());
    }

    // the first node is our root
    m_rootNode = internal.front();

    // for each internal node (including the root) assign one or two internal
    // node as its child[ren]
    vector<PhyloTreeNode*>::iterator it;
    for (it = internal.begin(); it != internal.end()-1; ++it) {
        if ((*(it+1))->isRoot()) {
            (*it)->addChild(*(it+1), rand()/(double)RAND_MAX);
        }

        // add one more with a probability of 0.5
        if (rand()%2==0 && it != internal.end()-2 && (*(it+2))->isRoot()) {
            (*it)->addChild(*(it+2), rand()/(double)RAND_MAX);
        }
    }

    // last, add the leaf nodes
    vector<PhyloTreeNode*>::iterator lt = leaves.begin();
    for (it = internal.begin(); it != internal.end(); ++it) {
        while ((*it)->numChildren() < 2) {
            (*it)->addChild(*(lt++), rand()/(double)RAND_MAX);
        }
    }

    // all leaves should now have a parent
    assert(lt == leaves.end());
}
