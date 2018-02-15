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
from diffsel_script_utils import param, data
import diffsel_script_utils as utils
from subprocess import Popen, PIPE, call
from os.path import isfile
from os import remove

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

def STEP(s):
    print(utils.step(s))

def INFO(s):
    print("-- "+s)

def ERROR(s):
    print("\n\t["+utils.boldred("Error")+"] "+s)

def BAD(s):
    print(utils.bad(s))

def GOOD(s):
    print(utils.good(s))

STEP("Parsing command line arguments")

args = docopt(__doc__)
for arg in args:
    if arg != "--help":
        INFO(strip(arg)+" set to "+param(to_str(args[arg])))

path = args["--tracecomp-path"][0]
chain1 = args["<trace1>"]
chain2 = args["<trace2>"]

STEP("Running tracecomp")

iterations = min(sum(1 for line in open(chain1+".trace")), sum(1 for line in open(chain2+".trace"))) - 1
INFO("Traces contain "+data(iterations)+" iterations")

burnin = int(int(args["--burnin"][0]) * iterations / 100.)
INFO("Burnin set to "+data(burnin)+" iterations")

INFO("Removing old tracecomp output file if present")
if isfile("./tmp.tracecomp.contdiff"):
    remove("./tmp.tracecomp.contdiff")

INFO("Running tracecomp with burnin : "+data(burnin))
process = Popen([path+"/tracecomp -x "+str(burnin)+" -o tmp.tracecomp "+chain1+" "+chain2], shell=True, stdout=PIPE, stderr=PIPE)
result = process.wait()

if result == 0 and isfile("tmp.tracecomp.contdiff"):
    INFO("Tracecomp executed; results stored in "+data("tmp.tracecomp.contdiff"))

    tracecomp_data = pd.read_csv("./tmp.tracecomp.contdiff", sep='\t')

    mineff = min(tracecomp_data.effsize)
    maxrel = max(tracecomp_data.rel_diff)
    INFO("Minimum effsize is "+data(mineff)+" and maximum rel_diff is "+data(maxrel))

    if mineff < 50 or maxrel > 0.3:
        BAD("The chains do not seem to have converged!")
    else:
        GOOD("The chains seem to have converged!")
        convergence = True
else:
    ERROR("Something went wrong with tracecomp!\n\tOutput of command was:\n")
    print(process.communicate()[1].decode("ascii"))


