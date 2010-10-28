#include "EvolutionModel.h"

#include <cassert>
#include <cmath>

using std::exp;

// helper function: returns true if given base is a purine
bool isPuri(char b)
{
    assert(b == 'A' || b == 'C' || b == 'G' || b == 'T');
    return (b == 'A' || b == 'G');
}

// helper function: returns true if given base is a pyrimidine
bool isPyri(char b)
{
    assert(b == 'A' || b == 'C' || b == 'G' || b == 'T');
    return (b == 'C' || b == 'T');
}

double Kimura::P(char a, char b, double t) const
{
    assert(isPuri(a) || isPyri(a));
    assert(isPuri(b) || isPyri(b));

    if (((isPuri(a) and isPuri(b)) or (isPyri(a) and isPyri(b))) and a != b) {
        // transition
        return 0.25-exp(-(2.*R+1.)/(R+1.)*t)/2.0+exp(-2.0/(R+1)*t)/4.0;
    }
    else if ((isPuri(a) and isPyri(b)) or (isPyri(a) and isPuri(b))) {
        // transversion
        return (0.5-exp(-2.0/(R+1)*t)/2.0)/2.0;
    } else {
        // no change of base
        return 0.25+exp(-(2.0*t)/(R+1))/4.0+exp(((-2.0*R-1)*t)/(R+1))/2.0;
    }
}
