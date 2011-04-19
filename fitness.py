#!/share/apps/modulessoftware/python/python-2.7/bin/python -u

from utils.genome2newick import *
import sys
import re
import subprocess
import tempfile
import os
import hashlib
import time

if len(sys.argv) < 2:
    sys.stderr.write("Usage: "+sys.argv[0]+" [phy file] <--no-lengths>\n");
    sys.exit(1)

lengths = not (len(sys.argv) == 3 and sys.argv[2] == "--no-lengths")

line = sys.stdin.readline()
while line:
    newick_data = newick(line, lengths)
    tree_file = '/tmp/'+hashlib.sha1(newick_data).hexdigest()
    tree_file = tempfile.NamedTemporaryFile(suffix='.nwk')

    tree_file.write(newick_data)
    tree_file.flush()

    args = ["phyml", "-i", sys.argv[1], "-m", "K80", "-f",
            "0.25,0.25,0.25,0.25", "-t", "10.0", "-u", tree_file.name, "-o",
            ("n" if lengths else "l")]

    start = time.time()
    output,_ = subprocess.Popen(args, stdout=subprocess.PIPE).communicate()

    m = re.search('Log likelihood of the current tree:\s+(-?\d+\.\d+)', output)
    logLH = m.group(1)

    sys.stderr.write(str(newick_data)+': '+logLH+'('+str(time.time()-start)+' s)\n')

    print logLH

    line = sys.stdin.readline()

sys.stderr.write("done\n")
