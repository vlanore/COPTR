"""EARLIEST_CONVERGENCE.PY
A script that, given two MCMC traces, tries to find the earliest iteration such that the
chains are considered convergent by tracecomp.

Usage:
  earliest_convergence.py [options...] <trace1> <trace2>

Options:
  -h, --help                           Show this message.
  -s <size>, --step-size <size>        Set the step size for the verification [default: 50].
  -t <path>, --tracecomp-path <path>   Set the path to tracecomp [default: ./].
  -b <burnin>, --burnin <burnin>       Set the burnin (as a percentage) [default: 20].
"""

import pandas as pd
from docopt import docopt
from diffsel_script_utils import *
from subprocess import Popen, PIPE, call
from os.path import isfile

def strip(s):
    if s[0] == '-' or s[0] == '<':
        return strip(s[1:])
    elif s[-1:] == '>':
        return s[:-1]
    else:
        return s

def to_str(thing):
    if isinstance(thing, int):
        return str(thing)
    elif not isinstance(thing, str):
        return str(thing[0])
    else:
        return str(thing)

print(step("Parsing command line arguments"))

args = docopt(__doc__)
for arg in args:
    if arg != "--help":
        print("-- "+strip(arg)+" set to "+param(to_str(args[arg])))

burnin = int(args["--burnin"][0])
path = args["--tracecomp-path"][0]
chain1 = args["<trace1>"]
chain2 = args["<trace2>"]

print(step("Running tracecomp"))

print("-- Running tracecomp with burnin : "+data(burnin))
process = Popen([path+"/tracecomp -x "+str(burnin)+" -o tmp.tracecomp "+chain1+" "+chain2], shell=True, stdout=PIPE, stderr=PIPE)
result = process.wait()

if result == 0 and isfile("tmp.tracecomp.contdiff"):
    print("-- Tracecomp executed; results stored in "+data("tmp.tracecomp.contdiff"))

    tracecomp_data = pd.read_csv("./tmp.tracecomp.contdiff", sep='\t')

    mineff = min(tracecomp_data.effsize)
    maxrel = max(tracecomp_data.rel_diff)
    print("-- Minimum effsize is "+data(mineff)+" and maximum rel_diff is "+data(maxrel))

    if mineff < 50 or maxrel > 0.3:
        print(bad("The chains do not seem to have converged!"))
    else:
        print(good("The chains seem to have converged!"))
        convergence = True
else:
    print("-- "+boldred("Warning")+" Something went wrong with tracecomp!\n-- Output of command was:\n")
    print(process.communicate()[1].decode("ascii"))


