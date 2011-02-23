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

/** concatenates a vector of strings.
 */
string detokenize(const vector<string>& tokens, string sep = string("\t"));

/** converts a string to a double.
 * \exception BadConversion String cannot be converted to float.
 */
double convertToDouble(string const& s);

#endif
