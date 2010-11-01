#ifndef EASYSTEM_H_INCLUDED
#define EASYSTEM_H_INCLUDED

#include <functional>
#include <vector>
#include <ostream>
#include <istream>
#include <string>
#include <cassert>
#include <sstream>

using std::vector;
using std::ostream;
using std::istream;
using std::string;
using std::binary_function;

// helper class for the runGenerations method: returns true after the given number
// of generations have been generated
template<typename T>
class Generations : public binary_function<vector<T>, int, bool>
{
    Generations(int n) : m_noGenerations(n) {}
    bool operator()(vector<T> _, int genNo) { return genNo >= m_noGenerations; }
private:
    int m_noGenerations;
};

template<typename T>
class EASystem {
public:
    EASystem() {}
    virtual ~EASystem() {}

    // the following two methods are made to facilitate the use of simdist
    void exportGenomes(ostream& out) const // write genomes to stdout
    {
        typename vector<T>::const_iterator it;
        for (it = m_population.begin(); it != m_population.end(); ++it) {
            out << *it;
        }
    }
    void readFitnessValues(istream& in) // read fitness values from stdin
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

    vector<T> getNBest(unsigned int n) const; // get the n best individuals

    void runUntil(binary_function<vector<T>, int, bool>& stoppingCriterion); // run until the given predicate, given fitness and gen. no, returns true
    void runGenerations(unsigned int noGenerations) // utility method: run for given number of generations
    {
        runUntil(Generations<T>(noGenerations));
    }

private:
    vector<T> m_population;
    vector<double> m_fitnessValues;
};

#endif
