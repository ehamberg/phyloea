#!/usr/bin/env python

import sys

def newick(line, branchLengths = True):
    lengths = []
    arity = [2]
    index = -1

    newick = ""

    for token in line.strip().split('\t'):
        index += 1

        if branchLengths and index % 2 == 1:
            lengths.append(token)
            continue

        arity[-1] -= 1
        if token == 'h':
            newick += '('
            arity.append(2)
        else:
            newick += token
            if branchLengths:
                newick += ':'
                newick += lengths.pop()
            if arity[-1] == 1:
                newick += ','

        while arity[-1] == 0 and (not branchLengths or len(lengths) > 0):
            newick += ')'
            if branchLengths:
                newick += ':' +  lengths.pop()
            del arity[-1]
            if arity[-1] == 1 and len(arity) > 1:
                newick += ','

    newick += ");"

    return newick

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--no-lengths":
        print(newick(sys.stdin.readline(), False))
    else:
        print(newick(sys.stdin.readline()))
