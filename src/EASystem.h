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

class lt {
    const vector<double>& _x;
public:
    lt(const vector<double>& x ) : _x(x) {}
    bool operator()( double j, double k ) const { return _x[j] > _x[k]; }
};

template<typename T>
class EASystem {
public:
    EASystem(MutationOp<T>* m, RecombOp<T>* r, SelectionOp<T> *s);
    virtual ~EASystem() { delete m_mutOp; delete m_recOp; delete m_selectionOp; }

    void exportGenomes(ostream& out) const; // write genomes to stdout
    void readFitnessValues(istream& in); // read fitness values from stdin
    vector<T> getNBest(unsigned int n) const;
    void runUntil(Generations<vector<T> > stoppingCriterion); // run until the given predicate, given fitness and gen. no, returns true
    void runGenerations(unsigned int noGenerations); // utility method: run for given number of generations
    void setPopulation(vector<T> pop);
    double averageFitness() const;
    double maxFitness() const;
    void setElitism(int e) { m_elitism = e; }
    void setLogStream(std::ostream* s) { m_logStream = s; }

private:
    vector<T> m_population;
    vector<double> m_fitnessValues;
    int m_elitism;
    int m_generationNumber;
    MutationOp<T>* m_mutOp;
    RecombOp<T>* m_recOp;
    SelectionOp<T>* m_selectionOp;
    std::ostream* m_logStream;
};

#endif
