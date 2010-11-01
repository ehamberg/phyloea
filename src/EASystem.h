#ifndef EASYSTEM_H_INCLUDED
#define EASYSTEM_H_INCLUDED

#include <functional>
#include <vector>
#include <ostream>
#include <istream>
#include <string>
#include <cassert>
#include <sstream>
#include <iostream>

using std::vector;
using std::ostream;
using std::istream;
using std::string;
using std::binary_function;

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
            out << *it << '\t';
        }
        out << '\n';
    }

    void readFitnessValues(istream& in) // read fitness values from stdin
    {
        string values;

        getline(in, values);

        std::cout << values << '\n';

        std::istringstream stm(values);

        double fitness;
        unsigned int i = 0;
        while (stm >> fitness) {
            std::cout << "***\n";
            m_fitnessValues[i++] = fitness;
        }

        assert(i == m_population.size());
        assert(m_population.size() == m_fitnessValues.size());
    }

    vector<T> getNBest(unsigned int n) const // get the n best individuals
    {
        // TODO:
        // Create a vector i = {0,1,2,3,4,5} of sort indexes (std::generate should
        // com in handy). Use the three argument sort, where the third argument
        // references your vector x, to compare the elements of i based not on
        // their own values but on the values they index in x:
        //
        // class lt {
        // vector<int&_x;
        // public:
        // lt( vector<int& x ) : _x(x) {}
        // bool comp( int j, int k ) const { return _x[j] < _x[k]; }
        // };
        // ....
        // sort( i.begin(), i.end(), lt(x) );
            }

    // FIXME
    void runUntil(Generations<vector<T> > stoppingCriterion) // run until the given predicate, given fitness and gen. no, returns true
    {
        while (!stoppingCriterion(m_population, m_generationNumber++)) {
            std::cout << "Generation " << m_generationNumber << '\n';
            exportGenomes(std::cout);
            readFitnessValues(std::cin);
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
        std::cout << m_population.size() << "lmao\n";
        std::cout << pop.size() << "lmao\n";
        m_fitnessValues.resize(m_population.size());
    }

private:
    vector<T> m_population;
    vector<double> m_fitnessValues;
    int m_generationNumber;
};

#endif
