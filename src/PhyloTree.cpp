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
}

PhyloTreeNode::PhyloTreeNode(string name, string states)
{
    this->name = name;
    this->states = states;
    this->nStates = states.size();
}

PhyloTreeNode::~PhyloTreeNode()
{
    //cout << "bye from " << name << '\n';
    //cout << "\tdeleting " << children.size() << " children..." << '\n';
    // delete all children
    vector<pair<PhyloTreeNode*, double> >::iterator it;
    for (it = children.begin(); it != children.end(); it++) {
        delete (*it).first;
    }
}

void PhyloTreeNode::addChild(PhyloTreeNode* child, double distance) {
    if (nStates == 0) {
        nStates = child->noStates();
    } else {
        assert(nStates == child->noStates());
    }
    children.push_back(pair<PhyloTreeNode*, double>(child, distance));
}

int PhyloTreeNode::height() const
{
    if (children.empty()) {
        return 1;
    }

    int maxH = 0;
    vector<pair<PhyloTreeNode*, double> >::const_iterator it;
    for (it = children.begin(); it != children.end(); it++) {
        maxH = MAX(maxH, 1+it->first->height());
    }

    return maxH;
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
    case 'G':
        likelihoods[1] = 1.0;
        break;
    case 'C':
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

    if (children.empty()) { // leaf node
        for (int i = 0; i < noStates(); i++) {
            likelihoods.push_back(leafLikelihood(states[i]));
        }

        return likelihoods;
    }

    double xs, ys;
    char nucleotides[4] = {'A', 'G', 'C', 'T' };

    for (int i = 0; i < noStates(); i++) {
        vector<double> l;
        for (int s = 0; s < 4; s++) {
            xs = ys = 0.0;
            for (int x = 0; x < 4; x++) {
                double p1 = eModel->P(nucleotides[s], nucleotides[x], children[0].second);
                double p2 = eModel->P(nucleotides[s], nucleotides[x], children[1].second);
                double l1 = children[0].first->likelihood(eModel)[i][x];
                double l2 = children[1].first->likelihood(eModel)[i][x];
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

    vector<pair<PhyloTreeNode*, double> >::const_iterator it;
    for (it = children.begin(); it != children.end(); it++) {
        out << "\t\"" << name << "\" -> \"" << it->first->getName()
            << "\" [label=\"" << it->second << "\"];\n";

        out << it->first->links();
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
