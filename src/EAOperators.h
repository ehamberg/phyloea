#ifndef EAOPERATORS_H_INCLUDED
#define EAOPERATORS_H_INCLUDED

#include <vector>

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
    SelectionOp();
    virtual ~SelectionOp();
    virtual T select(const std::vector<T>& pop, const std::vector<double> fitness) = 0;
};

// some pretty standard selection schemes
template<typename T>
class RouletteWheelSelection : public SelectionOp<T> {
public:
    RouletteWheelSelection();
    virtual ~RouletteWheelSelection();
    T select(const std::vector<T>& pop, const std::vector<double> fitness);
};

#endif
