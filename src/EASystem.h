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

// helper class for the runGenerations method: returns true after the given number
// of generations have been generated
template<typename T>
class Generations : public binary_function<T, int, bool> {
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
    EASystem(MutationOp<T>* m, RecombOp<T>* r, SelectionOp<T> *s)
        : m_mutOp(m), m_recOp(r), m_selectionOp(s) { m_elitism = 0; }
    virtual ~EASystem() { delete m_mutOp; delete m_recOp; delete m_selectionOp; }

    // the following two methods are made to facilitate the use of simdist
    void exportGenomes(ostream& out) const // write genomes to stdout
    {
        typename vector<T>::const_iterator it;
        for (it = m_population.begin(); it != m_population.end(); ++it) {
            out << *it << '\n';
        }
    }

    void readFitnessValues(istream& in) // read fitness values from stdin
    {
        string values;
        string temp;

        for (unsigned int i = 0; i < m_population.size(); ++i) {
            in >> temp;
            values += temp + '\n';
        }

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
            exportGenomes(std::cout);
            readFitnessValues(std::cin);

            vector<T> newGeneration;

            if (m_elitism > 0) {
                //std::cerr << "elitism: " << m_elitism << '\n';
                vector<T> elite = getNBest(m_elitism);
                newGeneration = elite;
            }

            typename vector<T>::const_iterator it;
            while (newGeneration.size() < m_population.size()) {
                T parent1 = m_selectionOp->select(m_population, m_fitnessValues);
                T parent2 = m_selectionOp->select(m_population, m_fitnessValues);

                // generate offspring ...
                vector<T> offspring = m_recOp->produceOffspring(parent1, parent2);

                // ... and add these to the new generation
                for (it = offspring.begin(); it != offspring.end(); ++it) {
                    newGeneration.push_back(*it);
                }
            }

            m_mutOp->mutate(newGeneration);
            m_population = newGeneration;

            std::cerr << "Generation " << m_generationNumber << "\tavg:\t" <<
                averageFitness() << "\tmax:\t" << highestFitness() << '\n';
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

    double averageFitness()
    {
        double sum = 0;
        for (cdit = m_fitnessValues.begin(); cdit != m_fitnessValues.end(); ++cdit) {
            sum += *cdit;
        }

        return sum/m_fitnessValues.size();
    }

    double highestFitness()
    {
        double highest = 0;
        for (cdit = m_fitnessValues.begin(); cdit != m_fitnessValues.end(); ++cdit) {
            if (*cdit > highest) {
                highest = *cdit;
            }
        }

        return highest;
    }

    void setElitism(int e) { m_elitism = e; }

private:
    vector<T> m_population;
    vector<double> m_fitnessValues;
    int m_elitism;
    int m_generationNumber;
    MutationOp<T>* m_mutOp;
    RecombOp<T>* m_recOp;
    SelectionOp<T>* m_selectionOp;
    vector<double>::const_iterator cdit;
};

#endif
