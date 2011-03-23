#!/usr/bin/env python

import sys

def newick(line):
    lengths = []
    arity = [2]
    index = -1

    newick = ""

    for token in line.split('\t'):
        index += 1

        if index % 2 == 1:
            lengths.append(token)
            continue

        arity[-1] -= 1
        if token == 'h':
            newick += '('
            arity.append(2)
        else:
            newick += token
            newick += ':'
            newick += lengths.pop()
            if arity[-1] == 1:
                newick += ','

        while arity[-1] == 0 and len(lengths) > 0:
            newick += '):'
            newick += lengths.pop()
            del arity[-1]
            if arity[-1] == 1 and len(arity) > 1:
                newick += ','

    newick += ");"

    return newick

if __name__ == "__main__":
    print(newick(sys.stdin.readline()))
