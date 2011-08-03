#include "AMaxOperators.h"
#include "utils.h"

#include <cstdlib>
#include <iostream>

void MutateString::mutate(std::vector<std::string>& genomes)
{
    for (auto it = genomes.begin(); it != genomes.end(); ++it) {
        for (unsigned int i = 0; i < it->length(); i++) {
            if (randZeroToOne() <= m_mutationRate) {
                (*it)[i] = (char)(rand()%2+(int)'0');
            }
        }
    }
}

std::vector<std::string> RecombineString::produceOffspring(const std::string& p1, const std::string& p2)
{
    std::vector<std::string> children;
    if (randZeroToOne() >= m_recombProb) {
        int p = rand()%p1.length(); // crossover point
        std::string child1 = p1.substr(0, p)+p2.substr(p);
        std::string child2 = p2.substr(0, p)+p1.substr(p);

        children.push_back(child1);
        children.push_back(child2);
    } else {
        children.push_back(p1);
        children.push_back(p2);
    }

    return children;
}
