#!/usr/bin/python -u

from utils.genome2newick import *
import sys
import re
import subprocess
import tempfile
import os

if len(sys.argv) < 2:
    sys.stderr.write("Usage: "+sys.argv[0]+" [phy file]\n");
    sys.exit(1)

line = sys.stdin.readline()
while line:
    newick_data = newick(line)
    tree_file = tempfile.NamedTemporaryFile(suffix='.nwk')
    #print newick_data
    tree_file.write(newick_data)
    tree_file.flush()

    #sys.stderr.write(tree_file.name+'\n')
    #print "----"
    #print(open(tree_file.name,"rb").read())
    #print "----"

    output,x = subprocess.Popen(["phyml", "-i", sys.argv[1], "-b", "4",
        "-m", "K80", "-f", "0.25,0.25,0.25,0.25", "-t", "10.0", "-u",
        tree_file.name, "-o", "n"], stdout=subprocess.PIPE).communicate()
    #sys.stderr.write(output+'\n')

    #print ["phyml", "-i", "tests/aligned.phy", "-b", "4", "-m", "K80", "-f",
    #       "0.25,0.25,0.25,0.25", "-t", "10.0", "-u", tree_file.name, "-o", "n"]

    m = re.search('Log likelihood of the current tree:\s+(-?\d+\.\d+)', output)
    logLH = m.group(1)

    print logLH

    line = sys.stdin.readline()

sys.stderr.write("done\n")
