#include "utils.h"

#include <sstream>

using std::istringstream;

vector<string> tokenize(const string& s)
{
    istringstream iss(s);
    vector<string> tokens;
    string temp;

    while (iss >> temp) tokens.push_back(temp);

    return tokens;
}

string detokenize(const vector<string>& tokens, string sep)
{
    string s;

    for (auto it = tokens.cbegin(); it != tokens.cend(); ++it) {
        s += *it;
        if (it != tokens.end()-1) {
            s += sep;
        }
    }

    return s;
}

double convertToDouble(string const& s)
{
    std::istringstream i(s);
    double x;
    if (!(i >> x))
        throw BadConversion("convertToDouble(\"" + s + "\")");
    return x;
}

string convertToString(double x)
{
    std::ostringstream o;
    if (!(o << x))
        throw BadConversion("stringify(double)");
    return o.str();
}

string convertToString(int x)
{
    std::ostringstream o;
    if (!(o << x))
        throw BadConversion("stringify(double)");
    return o.str();
}
