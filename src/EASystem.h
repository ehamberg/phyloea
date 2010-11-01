#ifndef EASYSTEM_H_INCLUDED
#define EASYSTEM_H_INCLUDED

#include <functional>
#include <vector>
#include <ostream>
#include <istream>
#include <string>

using std::vector;
using std::ostream;
using std::istream;
using std::string;
using std::binary_function;

template<typename T>
class EASystem {
public:
    EASystem();
    virtual ~EASystem();

    // the following two methods are made to facilitate the use of simdist
    void exportGenomes(ostream& out) const; // write genomes to stdout
    void readFitnessValues(istream& in); // read fitness values from stdin

    vector<T> getNBest(unsigned int n) const; // get the n best individuals

    void runUntil(binary_function<vector<T>, int, bool>& stoppingCriterion); // run until the given predicate, given fitness and gen. no, returns true
    void runGenerations(unsigned int noGenerations); // utility method: run for given number of generations

private:
    vector<T> m_population;
    vector<double> m_fitnessValues;
};

#endif
