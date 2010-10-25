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

    this->name = "node " + out.str();
}

PhyloTreeNode::PhyloTreeNode(string name, vector<char> states)
{
    this->name = name;
    this->states = states;
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
    children.push_back(pair<PhyloTreeNode*, double>(child, distance));
}

int PhyloTreeNode::height() const
{
    if (children.size() == 0) {
        return 1;
    }

    int maxH = 0;
    for (unsigned int i = 0; i < children.size(); i++) {
        maxH = MAX(maxH, 1+children[i].first->height());
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

string PhyloTreeNode::links() const
{
    stringstream out;

    for (unsigned int i = 0; i < children.size(); i++) {
        out << "\t\"" << name << "\" -> \"" << children[i].first->getName()
            << "\" [label=\"" << children[i].second << "\"];\n";

        out << children[i].first->links();
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
