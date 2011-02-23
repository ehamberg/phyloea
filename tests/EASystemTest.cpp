#include "gtest/gtest.h"

#include <vector>
#include <string>
#include <cmath>
#include <sstream>

#include "EASystem.h"
#include "EAOperators.h"

class MutateString : public MutationOp<std::string> {
public:
    MutateString(double mutRate) : m_mutationRate(mutRate) {}
    ~MutateString() {}
    void mutate(std::vector<std::string>& genomes)
    {
        std::vector<std::string>::iterator it;
        for (it = genomes.begin(); it != genomes.end(); ++it) {
            for (unsigned int i = 0; i < it->length(); i++) {
                if ((double)rand()/(double)RAND_MAX <= m_mutationRate) {
                    (*it)[i] = (char)(rand()%2+(int)'0');
                }
            }
        }
    }

private:
    double m_mutationRate;
};

class RecombineString : public RecombOp<std::string> {
public:
    RecombineString(double recombProb) : m_recombProb(recombProb) {}
    ~RecombineString() {}
    std::vector<std::string> produceOffspring(const std::string& p1, const std::string& p2)
    {
        std::vector<std::string> children;
        if ((double)rand()/(double)RAND_MAX >= m_recombProb) {
            int p = rand()%p1.length(); // crossover point
            std::string child1 = p1.substr(0, p)+p2.substr(p);
            std::string child2 = p2.substr(0, p)+p1.substr(p);

            children.push_back(child1);
            children.push_back(child2);
        } else {
            children.push_back(p1);
            children.push_back(p2);
        }

        return children;
    }
private:
    double m_recombProb;
};

class CountOnes : public FitnessFunc<std::string> {
    std::vector<double> fitness(const std::vector<std::string>& genomes)
    {
        std::vector<double> fitness;
        std::vector<std::string>::const_iterator it;
        for (it = genomes.begin(); it != genomes.end(); ++it) {
            double f = 0;
            for (unsigned int i = 0; i < it->length(); i++) {
                if (it->at(i) == '1') {
                    f++;
                }
            }
            fitness.push_back(f);
        }
        return fitness;
    }
};

class EASystemTest : public ::testing::Test {
protected:
    EASystemTest() {
        ea = new EASystem<std::string>(new MutateString(0.1),
                new RecombineString(0.7), new RouletteWheelSelection<string>,
                new CountOnes);

        std::vector<std::string> randomPop;
        for (int i = 0; i < 10; i++) {
            string s;
            for (int j = 0; j < 20; j++) {
                s += (char)(rand()%2+(int)'0');
            }
            randomPop.push_back(s);
        }

        ea->setPopulation(randomPop);
    }

    virtual ~EASystemTest() {
        delete ea;
    }

    EASystem<std::string>* ea;
};

bool monotonicRising(const std::vector<double>& v)
{
    for (unsigned int i = 0; i < v.size()-1; i++) {
        if (v[i] > v[i+1]) {
            return false;
        }
    }

    return true;
}

// Make sure max fitness never decreases when using elitism
TEST_F(EASystemTest, MonotonicWithElitism) {
    ea->setElitism(4);
    std::ostringstream oss;

    ea->setLogStream(&oss);
    ea->runGenerations(50);

    std::istringstream iss(oss.str());

    double s;
    int i = 0;
    std::vector<double> fitnessHist;
    while (iss >> s) {
        if (i++%4 == 2) {
            fitnessHist.push_back(s);
        }
    }

    ASSERT_TRUE(monotonicRising(fitnessHist));
}
