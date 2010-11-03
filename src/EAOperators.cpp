#include "EAOperators.h"

#include <cassert>
#include <cstdlib>

template<typename T>
T RouletteWheelSelection<T>::select(const std::vector<T>& pop, const std::vector<double> fitness)
{
    assert(pop.size() > 0);

    double fitnessSum = 0;
    std::vector<double>::const_iterator it;

    for (it = fitness.begin(); it != fitness.end(); it++) {
        fitnessSum += *it;
    }

    // get a random number between 0 and fitness sum ...
    double r = (double)rand()/(double)RAND_MAX*fitnessSum;

    // ... then run through fitness values until sum >= fitnessSum and pick
    // the individual in that position
    int i;
    for (i = 0; fitnessSum >= fitness[i]; i++) {
        fitnessSum -= fitness[i];
    }

    return pop[i];
}
