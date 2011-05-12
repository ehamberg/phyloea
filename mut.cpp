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

    for (int i = 0; i < 99000; i++) {
        vector<string> genomes;
        genomes.push_back(genome);
        MutateTree mt(20, 1.0);
        mt.mutate(genomes);
        cout << genomes.at(0) << endl;
    }

    return 0;
}
