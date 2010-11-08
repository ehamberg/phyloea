#include "PhyloTree.h"

#include <sstream>
#include <iostream>
#include <cassert>
#include <cmath>

using std::vector;
using std::pair;
using std::string;
using std::stringstream;
using std::cout;

#define MAX(x,y) (x>y?x:y)

int PhyloTreeNode::count = 0;

PhyloTreeNode::PhyloTreeNode()
{
    // generate a name for this node of the form “node #”
    stringstream out;
    out << PhyloTreeNode::count++;
    out.str();

    m_name = "node " + out.str();
    m_nStates = 0;

    m_left = m_right = NULL;
}

PhyloTreeNode::PhyloTreeNode(string name, string states)
{
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

    double xs, ys;
    char nucleotides[4] = {'A', 'C', 'G', 'T' };

    for (unsigned int i = 0; i < noStates(); i++) {
        vector<double> l;
        for (int s = 0; s < 4; s++) {
            xs = ys = 0.0;
            for (int x = 0; x < 4; x++) {
                double p1 = eModel->P(nucleotides[s], nucleotides[x], m_leftDist);
                double p2 = eModel->P(nucleotides[s], nucleotides[x], m_rightDist);
                double l1 = m_left->likelihood(eModel)[i][x];
                double l2 = m_right->likelihood(eModel)[i][x];
                xs += p1*l1;
                ys += p2*l2;
            }
            l.push_back(xs*ys);
        }

        m_likelihoods.push_back(l);
    }

    return m_likelihoods;
}

// print all links from this node to child nodes in graphviz format
string PhyloTreeNode::links() const
{
    stringstream out;

    if (m_left) {
        out << "\t\"" << m_name << "\" -> \"" << m_left->getName()
            << "\" [label=\"" << m_leftDist << "\"];\n";
        out << m_left->links();
    }

    if (m_right) {
        out << "\t\"" << m_name << "\" -> \"" << m_right->getName()
            << "\" [label=\"" << m_leftDist << "\"];\n";
        out << m_right->links();
    }

    return out.str();
}

ostream& operator<<(ostream& out, const PhyloTree& t)
{
    out << t.m_rootNode->links();

    return out;
}

string PhyloTree::dot() const
{
    return "digraph {\n" + m_rootNode->links() + "}";
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
