#ifndef EASYSTEM_H_INCLUDED
#define EASYSTEM_H_INCLUDED

#include <functional>
#include <algorithm>
#include <vector>
#include <ostream>
#include <istream>
#include <string>
#include <cassert>
#include <sstream>
#include <iostream>

#include "EAOperators.h"

using std::vector;
using std::ostream;
using std::istream;
using std::string;
using std::binary_function;

void printVal(string i) {
    std::cout << i << '\n';
}

// helper class for the runGenerations method: returns true after the given number
// of generations have been generated
template<typename T>
class Generations : public binary_function<T, int, bool>
{
public:
    Generations(int n) : m_noGenerations(n) {}
    bool operator()(T _, int genNo) { return genNo >= m_noGenerations; }
private:
    int m_noGenerations;
};

// class generator:
struct c_unique {
  int current;
  c_unique() {current=0;}
  int operator()() {return current++;}
} UniqueNumber;

class lt {
    const vector<double>& _x;
public:
    lt(const vector<double>& x ) : _x(x) {}
    bool operator()( double j, double k ) const { return _x[j] > _x[k]; }
};

template<typename T>
class EASystem {
public:
    EASystem(MutationOp<T>* m, RecombOp<T>* r) : m_mutOp(m), m_recOp(r) {}
    virtual ~EASystem() { delete m_mutOp; delete m_recOp; }

    // the following two methods are made to facilitate the use of simdist
    void exportGenomes(ostream& out) const // write genomes to stdout
    {
        typename vector<T>::const_iterator it;
        for (it = m_population.begin(); it != m_population.end(); ++it) {
            out << *it << '\t';
        }
        out << '\n';
    }

    void readFitnessValues(istream& in) // read fitness values from stdin
    {
        string values;

        getline(in, values);

        std::cout << "read fitvals: " << values << '\n';

        std::istringstream stm(values);

        double fitness;
        unsigned int i = 0;
        while (stm >> fitness) {
            m_fitnessValues[i++] = fitness;
        }

        assert(i == m_population.size());
        assert(m_population.size() == m_fitnessValues.size());
    }

    vector<T> getNBest(unsigned int n) const // get the n best individuals
    {
        assert(n <= m_population.size());
        vector<int> indices(m_population.size());
        vector<int>::iterator it;

        std::generate(indices.begin(), indices.end(), UniqueNumber);
        sort(indices.begin(), indices.end(), lt(m_fitnessValues));

        vector<T> ret(n);
        for (unsigned int i = 0; i < n; i++) {
            ret.push_back(m_population[indices[i]]);
        }

        return ret;
    }

    void runUntil(Generations<vector<T> > stoppingCriterion) // run until the given predicate, given fitness and gen. no, returns true
    {
        while (!stoppingCriterion(m_population, m_generationNumber++)) {
            std::cout << "Generation " << m_generationNumber << '\n';
            exportGenomes(std::cout);
            readFitnessValues(std::cin);

            vector<T> newGeneration(m_population.size());

            while (newGeneration.size() < m_population.size()) {
            }

            m_population = newGeneration;

            //for_each(foo.begin(), foo.end(), printVal);
        }
    }

    void runGenerations(unsigned int noGenerations) // utility method: run for given number of generations
    {
        runUntil(Generations<vector<T> >(noGenerations));
    }

    void setPopulation(vector<T> pop)
    {
        m_population = pop;
        m_generationNumber = 0;
        m_fitnessValues.resize(m_population.size());
    }

private:
    vector<T> m_population;
    vector<double> m_fitnessValues;
    int m_generationNumber;
    MutationOp<T>* m_mutOp;
    RecombOp<T>* m_recOp;
};

#endif
