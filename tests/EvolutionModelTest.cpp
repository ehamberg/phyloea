#include "gtest/gtest.h"
#include <cstdlib>
#include <cmath>

#include "EvolutionModel.h"

class EvolutionModelTest : public ::testing::Test {
protected:
    EvolutionModelTest() {
        srand(time(NULL));
        nucleotides[0] = 'A';
        nucleotides[1] = 'G';
        nucleotides[2] = 'C';
        nucleotides[3] = 'T';
    }

    virtual ~EvolutionModelTest() {
    }

    char nucleotides[4];
};

// test a simple tree for which the likelihood is known
TEST_F(EvolutionModelTest, Sums_to_one) {
    Kimura kimura((rand()%100)/10.0);
    double t = rand();
    int n = rand()%4;

    double sum = 0;
    for (int i = 0; i < 4; i++) {
        sum += kimura.P(nucleotides[n], nucleotides[i], t);
    }

    ASSERT_FLOAT_EQ(sum, 1.0);
}

TEST_F(EvolutionModelTest, P_less_than_one) {
    Kimura kimura((rand()%100)/10.0);

    for (int i = 0; i < 1000; i++) {
        double t = rand();
        int n1 = rand()%4;
        int n2 = rand()%4;
        double R = (rand()%100)/10.0;

        kimura.setR(R);

        ASSERT_LE(kimura.P(nucleotides[n1], nucleotides[n2], t), 1.0);
    }
}
