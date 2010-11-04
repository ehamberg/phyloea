#include "EASystem.h"

// generator for number sequence 0, 1, 2, …
struct c_unique {
  int current;
  c_unique() {current = 0;}
  int operator()() {return current++;}
} NextNum;

template<typename T>
void print(T s) {
    std::cerr << s << '\n';
}

template <typename T>
EASystem<T>::EASystem(MutationOp<T>* m, RecombOp<T>* r, SelectionOp<T> *s)
{
    m_elitism = 0;
    m_mutOp = m;
    m_recOp = r;
    m_selectionOp = s;
    m_logStream = NULL;
}

// the following two methods are made to facilitate the use of simdist
template <typename T>
void EASystem<T>::exportGenomes(ostream& out) const // write genomes to stdout
{
    typename vector<T>::const_iterator it;
    for (it = m_population.begin(); it != m_population.end(); ++it) {
        out << *it << '\n';
    }
}

template<typename T>
void EASystem<T>::readFitnessValues(istream& in) // read fitness values from stdin
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

template <typename T>
vector<T> EASystem<T>::getNBest(unsigned int n) const // get the n best individuals
{
    assert(n <= m_population.size());
    vector<int> indices(m_population.size());

    // generate a list of 0..n and sort according to the values in m_fitnessValues
    std::generate(indices.begin(), indices.end(), NextNum);
    sort(indices.begin(), indices.end(), lt(m_fitnessValues));

    vector<T> elite;
    elite.reserve(4);

    // push the n best individuals to the return vector
    for (unsigned int i = 0; i < n; i++) {
        elite.push_back(m_population[indices[i]]);
    }

    return elite;
}

// run until the given predicate, given fitness and gen. no, returns true
template <typename T>
void EASystem<T>::runUntil(Generations<vector<T> > stoppingCriterion)
{
    while (!stoppingCriterion(m_population, m_generationNumber++)) {
        exportGenomes(std::cout);
        readFitnessValues(std::cin);

        vector<T> newGeneration;

        typename vector<T>::const_iterator it;
        while (newGeneration.size()+m_elitism < m_population.size()) {
            T parent1 = m_selectionOp->select(m_population, m_fitnessValues);
            T parent2 = m_selectionOp->select(m_population, m_fitnessValues);

            // generate offspring ...
            vector<T> offspring = m_recOp->produceOffspring(parent1, parent2);

            // ... and add these to the new generation
            for (it = offspring.begin(); it != offspring.end(); ++it) {
                newGeneration.push_back(*it);
            }
        }

        m_mutOp->mutate(newGeneration); // FIXME: lmao, don't mutate the elite!

        // if we have elitism, add the best individuals to the next gen. directly
        if (m_elitism > 0) {
            vector<T> elite = getNBest(m_elitism);
            for (unsigned int i = 0; i < elite.size(); i++) {
                newGeneration.push_back(elite[i]);
            }
        }

        m_population = newGeneration;

        if (m_logStream) {
            *m_logStream << averageFitness() << '\t' << maxFitness() << '\n';
        }
        std::cerr << "Generation " << m_generationNumber << "\tavg:\t" <<
            averageFitness() << "\tmax:\t" << maxFitness() << '\n';
    }
}

// utility method: run for given number of generations
template <typename T>
void EASystem<T>::runGenerations(unsigned int noGenerations) 
{
    runUntil(Generations<vector<T> >(noGenerations));
}

template <typename T>
void EASystem<T>::setPopulation(vector<T> pop)
{
    m_population = pop;
    m_generationNumber = 0;
    m_fitnessValues.resize(m_population.size());
}

template <typename T>
double EASystem<T>::averageFitness() const
{
    double sum = 0;
    vector<double>::const_iterator cdit;
    for (cdit = m_fitnessValues.begin(); cdit != m_fitnessValues.end(); ++cdit) {
        sum += *cdit;
    }

    return sum/m_fitnessValues.size();
}

template <typename T>
double EASystem<T>::maxFitness() const // TODO : const
{
    double highest = 0;
    vector<double>::const_iterator cdit;
    for (cdit = m_fitnessValues.begin(); cdit != m_fitnessValues.end(); ++cdit) {
        if (*cdit > highest) {
            highest = *cdit;
        }
    }

    return highest;
}

// force methods for std::string
template EASystem<std::string>::EASystem(MutationOp<std::string>* m, RecombOp<std::string>* r, SelectionOp<std::string> *s);
template void EASystem<std::string>::exportGenomes(ostream& out) const;
template void EASystem<std::string>::setPopulation(vector<std::string> pop);
template void EASystem<std::string>::runGenerations(unsigned int noGenerations);
