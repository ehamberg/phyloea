#ifndef EAOPERATORS_H_INCLUDED
#define EAOPERATORS_H_INCLUDED

#include <vector>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <algorithm>

class lt;

template<typename T>
class MutationOp {
public:
    MutationOp() {}
    virtual ~MutationOp() {}
    virtual void mutate(std::vector<T>& genomes) = 0;
};

template<typename T>
class RecombOp {
public:
    RecombOp() {}
    virtual ~RecombOp() {}
    virtual std::vector<T> produceOffspring(const T& parent1, const T& parent2) = 0;
};

template<typename T>
class SelectionOp {
public:
    SelectionOp() {}
    virtual ~SelectionOp() {}
    virtual const T& select(const std::vector<T>& pop, const std::vector<double> fitness) = 0;
};

template <typename T>
class FitnessFunc {
public:
    FitnessFunc() {}
    virtual ~FitnessFunc() {}
    virtual std::vector<double> fitness(const std::vector<T>& genomes) = 0;
};


// some pretty standard selection schemes
template<typename T>
class RouletteWheelSelection : public SelectionOp<T> {
public:
    RouletteWheelSelection() {}
    virtual ~RouletteWheelSelection() {}
    const T& select(const std::vector<T>& pop, const std::vector<double> fitness)
    {
        assert(pop.size() > 0);

        double fitnessSum = 0;
        std::vector<double>::const_iterator it;

        for (it = fitness.begin(); it != fitness.end(); it++) {
            fitnessSum += *it;
        }

        // if total fitness is 0, return a random individual
        if (fitnessSum == 0.0) {
            return pop[rand()%pop.size()];
        }

        // get a random number between 0 and fitness sum ...
        double r = (double)rand()/(double)RAND_MAX*fitnessSum;

        // ... then run through fitness values until sum >= fitnessSum and pick
        // the individual in that position
        int i;
        for (i = 0; r >= fitness[i]; i++) {
            r -= fitness[i];
        }

        return pop[i];
    }
};

template<typename T>
class RankSelection : public SelectionOp<T> {
public:
    RankSelection() {}
    virtual ~RankSelection() {}
    const T& select(const std::vector<T>& pop, const std::vector<double> fitness)
    {
        assert(pop.size() > 0);

        std::vector<unsigned int> indices;
        unsigned int n = fitness.size();
        unsigned int rankSum =  n/2.0*(n+1);

        for (unsigned int i = 0; i < n; i++) {
            indices.push_back(i);
        }
        sort(indices.begin(), indices.end(), lt(fitness));

        std::vector<unsigned int> ranks(n);
        for (unsigned int i = 0; i < indices.size(); i++) {
            ranks[indices[i]] = indices.size()-i;
        }

        // get a random number between 0 and rank sum ...
        unsigned int r = rand()%rankSum;

        // ... then run through fitness values until sum >= rankSum and pick
        // the individual in that position
        int i = 0;
        for (; r >= ranks[i]; i++) {
            r -= ranks[i];
        }

        return pop[i];
    }
};

// for use with simdist: print
template<typename T>
class PipeFitnessFunc : public FitnessFunc<T> {
public:
    PipeFitnessFunc() {}
    virtual ~PipeFitnessFunc() {}
    virtual std::vector<double> fitness(const std::vector<T>& genomes)
    {
        std::vector<double> fitnessVals(genomes.size());

        typename std::vector<T>::const_iterator it;
        int ii = 1;
        for (it = genomes.begin(); it != genomes.end(); ++it) {
            std::cerr << ii++ << " requesting fitness for " << *it << '\n';
            std::cout << *it << '\n';
        }

        std::string values;
        std::string temp;

        for (unsigned int i = 0; i < genomes.size(); ++i) {
            std::cin >> temp;
            values += temp + '\n';
        }

        std::istringstream stm(values);

        double fitness;
        unsigned int i = 0;
        while (stm >> fitness) {
            fitnessVals[i++] = fitness;
        }

        assert(i == genomes.size());
        assert(genomes.size() == fitnessVals.size());

        return fitnessVals;
    }
};

#endif
