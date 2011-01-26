#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <vector>
#include <string>
#include <stdexcept>

using std::vector;
using std::string;
using std::runtime_error;

class BadConversion : public runtime_error {
public:
    BadConversion(string const& s) : runtime_error(s) {}
};

/** takes a string and returns a vector of its (whitespace-separated) sub
 * strings.
 */
vector<string> tokenize(const string& s);

double convertToDouble(string const& s);

string detokenize(const vector<string>& tokens, string sep = string("\t"));

#endif
