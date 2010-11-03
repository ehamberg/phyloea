#ifndef EAOPERATORS_H_INCLUDED
#define EAOPERATORS_H_INCLUDED

#include <vector>
#include <cassert>
#include <cstdlib>

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

#endif
