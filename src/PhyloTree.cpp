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

    name = "node " + out.str();
    nStates = 0;

    left = right = NULL;
}

PhyloTreeNode::PhyloTreeNode(string name, string states)
{
    this->name = name;
    this->states = states;
    this->nStates = states.size();

    left = right = NULL;
}

PhyloTreeNode::~PhyloTreeNode()
{
    // delete children
    delete left;
    delete right;
}

void PhyloTreeNode::addChild(PhyloTreeNode* child, double distance) {
    assert(!left or !right); // only two child nodes are allowed

    if (nStates == 0) {
        nStates = child->noStates();
    } else {
        assert(nStates == child->noStates());
    }

    if (!left) {
        left = child;
        leftDist = distance;
    } else if (!right) {
        right = child;
        rightDist = distance;
    }
}

int PhyloTreeNode::height() const
{
    if (!left and !right) {
        return 1;
    }

    if (!right) {
        return 1+left->height();
    }

    return 1+MAX(left->height(), right->height());
}

int PhyloTree::height() const
{
    return rootNode->height();
}

vector<double> PhyloTreeNode::leafLikelihood(char n) const {
    vector<double> likelihoods;
    for (int i = 0; i < 4; i++) {
        likelihoods.push_back(0.0);
    }

    switch (n) {
    case 'A':
        likelihoods[0] = 1.0;
        break;
    case 'C':
        likelihoods[1] = 1.0;
        break;
    case 'G':
        likelihoods[2] = 1.0;
        break;
    case 'T':
        likelihoods[3] = 1.0;
        break;
    }

    return likelihoods;
}

vector<vector<double> > PhyloTreeNode::likelihood(EvolutionModel* eModel)
{
    //  base case: we have already found likelihoods for the node
    if (!likelihoods.empty()) {
        return likelihoods;
    }

    if (!left and !right) { // leaf node
        for (unsigned int i = 0; i < noStates(); i++) {
            likelihoods.push_back(leafLikelihood(states[i]));
        }

        return likelihoods;
    }

    double xs, ys;
    char nucleotides[4] = {'A', 'C', 'G', 'T' };

    for (unsigned int i = 0; i < noStates(); i++) {
        vector<double> l;
        for (int s = 0; s < 4; s++) {
            xs = ys = 0.0;
            for (int x = 0; x < 4; x++) {
                double p1 = eModel->P(nucleotides[s], nucleotides[x], leftDist);
                double p2 = eModel->P(nucleotides[s], nucleotides[x], rightDist);
                double l1 = left->likelihood(eModel)[i][x];
                double l2 = right->likelihood(eModel)[i][x];
                xs += p1*l1;
                ys += p2*l2;
            }
            l.push_back(xs*ys);
        }

        likelihoods.push_back(l);
    }

    return likelihoods;
}

// print all links from this node to child nodes in graphviz format
string PhyloTreeNode::links() const
{
    stringstream out;

    if (left) {
        out << "\t\"" << name << "\" -> \"" << left->getName()
            << "\" [label=\"" << leftDist << "\"];\n";
        out << left->links();
    }

    if (right) {
        out << "\t\"" << name << "\" -> \"" << right->getName()
            << "\" [label=\"" << leftDist << "\"];\n";
        out << right->links();
    }

    return out.str();
}

ostream& operator<<(ostream& out, const PhyloTree& t)
{
    out << t.rootNode->links();

    return out;
}

string PhyloTree::dot() const
{
    return "digraph {\n" + rootNode->links() + "}";
}

double PhyloTree::logLikelihood()
{
    double ret = 0.0;

    vector<vector<double> > siteLikelihoods = rootNode->likelihood(evModel);

    vector<vector<double> >::const_iterator it;
    vector<double>::const_iterator jt;
    for (it = siteLikelihoods.begin(); it != siteLikelihoods.end(); it++) {

        double siteSum = 0.0;
        for (jt = (*it).begin(); jt != (*it).end(); jt++) {
            siteSum += (*jt);
        }

        ret += log10(siteSum/4.0);
    }

    return ret;
}

PhyloTree::~PhyloTree()
{
    delete evModel;
    delete rootNode; // FIXME
}
