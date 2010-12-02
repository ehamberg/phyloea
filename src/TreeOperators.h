#ifndef TREEOPERATORS_H_INCLUDED
#define TREEOPERATORS_H_INCLUDED

#include"EAOperators.h"

#include <vector>
#include <string>

class MutateTree : public MutationOp<std::string> {
public:
    MutateTree(unsigned int numNodes, double mutRate) : m_mutationRate(mutRate), m_numNodes(numNodes) {}
    ~MutateTree() {}
    void mutate(std::vector<std::string>& genome);
private:
    double m_mutationRate;
    unsigned int m_numNodes;
};

class RecombineTree : public RecombOp<std::string> {
public:
    RecombineTree(double recombProb) : m_recombProb(recombProb) {}
    ~RecombineTree() {}
    std::vector<std::string> produceOffspring(const std::string& p1, const std::string& p2);
private:
    double m_recombProb;
};

#endif
