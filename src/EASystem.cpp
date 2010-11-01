#include "EASystem.h"
#include <cassert>
#include <sstream>

template<typename T>
class Generations : public binary_function<vector<T>, int, bool>
{
    Generations(int n) : m_noGenerations(n) {}
    bool operator()(vector<T> _, int genNo) { return genNo >= m_noGenerations; }
private:
    int m_noGenerations;
};

template<typename T>
void EASystem<T>::exportGenomes(ostream& out) const
{
    typename vector<T>::const_iterator it;
    for (it = m_population.begin(); it != m_population.end(); ++it) {
        out << *it;
    }
}

template<typename T>
void EASystem<T>::readFitnessValues(istream& in)
{
    string values;

    in >> values;

    std::istringstream stm(values);

    double fitness;
    unsigned int i = 0;
    while (stm >> fitness) {
        m_fitnessValues[i++] = fitness;
    }

    assert(i == m_population.size()-1);
    assert(m_population.size() == m_fitnessValues.size());
}

template<typename T>
vector<T> EASystem<T>::getNBest(unsigned int n) const
{
}

template<typename T>
void EASystem<T>::runUntil(binary_function<vector<T>, int, bool>& stoppingCriterion)
{
}

template<typename T>
void EASystem<T>::runGenerations(unsigned int noGenerations)
{
    runUntil(Generations<T>(noGenerations));
}
