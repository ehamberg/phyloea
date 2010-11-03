#ifndef AMAXOPERATORS_H_INCLUDED
#define AMAXOPERATORS_H_INCLUDED

#include"EAOperators.h"

#include <vector>
#include <string>

class MutateString : public MutationOp<std::string> {
public:
    MutateString(double mutRate) : m_mutationRate(mutRate) {}
    ~MutateString() {}
    void mutate(std::vector<std::string>& genome);
private:
    double m_mutationRate;
};

class RecombineString : public RecombOp<std::string> {
public:
    RecombineString(double recombProb) : m_recombProb(recombProb) {}
    ~RecombineString() {}
    std::vector<std::string> produceOffspring(const std::string& p1, const std::string& p2);
private:
    double m_recombProb;
};

#endif
