#include "PhyloTree.h"
#include "EvolutionModel.h"
#include "Fasta.h"
#include "EASystem.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

#include "EAOperators.h"
#include "TreeOperators.h"

using namespace std;

int main(int argc, const char *argv[])
{
    string genome;
    getline(cin, genome);

    RecombineTree mt(1.0);
    for (int i = 0; i < 10; i++) {
        vector<string> offspring = mt.produceOffspring(genome, genome);
        cout << offspring.at(0) << endl;
    }

    return 0;
}
