# Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2017/11/27)
# Contributors:
# - Vincent Lanore <vincent.lanore@gmail.com>

# This software is a computer program whose purpose is to provide small tools and scripts related to phylogeny and bayesian
# inference.

# This software is governed by the CeCILL-B license under French law and abiding by the rules of distribution of free software.
# You can use, modify and/ or redistribute the software under the terms of the CeCILL-B license as circulated by CEA, CNRS and
# INRIA at the following URL "http://www.cecill.info".

# As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users
# are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive
# licensors have only limited liability.

# In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or
# reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated
# to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth
# computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements
# in conditions enabling the security of their systems and/or data to be ensured and, more generally, to use and operate it in
# the same conditions as regards security.

# The fact that you are presently reading this means that you have had knowledge of the CeCILL-B license and that you accept
# its terms.

from diffsel_script_utils import *
from os.path import isfile

#===================================================================================================
print(step("Parsing command line arguments"))

from argparse import ArgumentParser, FileType
parser = ArgumentParser(description='A non-regression test for diffsel that generates a toy alignment+tree, runs two diffsel chains and analyzes the results.')
parser.add_argument('-n', '--name', default="tmp_nonreg", help="the name of the run (used to name output files); default is 'tmp_nonreg'")
parser.add_argument('-u', '--until', type=int, default=10000, help="number of iterations to run diffsel for; default is 10k")
parser.add_argument('-d', '--depth', type=int, default=7, help="depth of the phylogenetic tree (which is binary and balanced); default is 7")
parser.add_argument('-b', '--burnin', type=float, default=0.2, help="burnin (as a float between 0 and 1); default is 0.2")
parser.add_argument('-a', '--analysis-only', dest="analysis", action='store_true', help="set this flag to run only the analysis part of the test\
 (you need to have diffsel traces already)")
parser.add_argument('-f', '--failures-allowed', dest="failure", action='store_true', help="set this flag to ignore (ie, not return an error) if the test fails")
args = parser.parse_args()

chainname = args.name
print("-- Chain name is "+param(chainname));
iterations = args.until
print("-- Number of iterations is "+param(iterations))
depth = args.depth
print("-- Depth of the tree is "+param(depth))
burninin = args.burnin
print("-- Burn-in is "+param(burninin * 100)+"%")
analysis = args.analysis
print("-- Analysis only: "+param(analysis))
allowfailure = args.failure
print("-- Failure allowed: "+param(analysis))

from subprocess import Popen, PIPE, call
if not analysis:
    #===================================================================================================
    print(step("Creating tree and alignment files"))

    print("-- Generating tree and alignment... ", end='')
    p = call(["python3 script/generate_toy "+str(depth)+" "+chainname], stdout=PIPE, shell=True)
    print("done")

    # print(p.stdout.readline())

    #===================================================================================================
    print(step("Launching diffsel chains"))

    print("-- Launching first chain")
    p2 = Popen(["time _build/diffsel -d "+chainname+".ali -t "+chainname+".tree -ncond 2 -x 1 "+str(iterations)+" "+chainname], stderr=PIPE, stdout=PIPE, shell=True)

    print("-- Launching second chain")
    p3 = Popen(["time _build/diffsel -d "+chainname+".ali -t "+chainname+".tree -ncond 2 -x 1 "+str(iterations)+" "+chainname+"2"], stderr=PIPE, stdout=PIPE, shell=True)

    print("-- Launching done; starting chain monitoring:")
    try:
        import progressbar
        bar = progressbar.ProgressBar()
        myrange = bar(range(iterations))
    except:
        myrange = range(iterations)
        print("-- Progressbar is not installed! Cannot display progress. Please wait for a while...")

    for i in myrange:
        p2.stdout.readline()
        p3.stdout.readline()

    r2 = p2.wait()
    r3 = p3.wait()
    if r2 == 0 and r3 == 0:
        print("-- Diffsel chains have ended succesfully!")
    else:
        print("-- "+boldred("Error")+": There was an error while running diffsel!\n-- Output of the first chain was:\n")
        print(p2.communicate()[1].decode("ascii"))
        exit(1)

#===================================================================================================
print(step("Chain convergence analysis"))
convergence = False

import os.path
if isfile(chainname+".trace") and isfile(chainname+"2.trace"):
    print("-- Found trace files: "+data(chainname+".trace")+" and "+data(chainname+"2.trace"))
else:
    print("-- "+boldred("Error")+": trace files do not exist! Expected files "+data(chainname+".trace")+" and "+data(chainname+"2.trace"))
    exit(1)

burnin = int(iterations * burninin)
print("-- Running tracecomp with burnin : "+data(burnin))
p4 = Popen(["_build/tracecomp -x "+str(burnin)+" -o "+chainname+" "+chainname+" "+chainname+"2"], shell=True, stdout=PIPE, stderr=PIPE)
r4 = p4.wait()

if r4 == 0 and isfile(chainname+".contdiff"):
    print("-- Tracecomp executed; results stored in "+data(chainname+".contdiff"))
    tcout = p4.communicate()[0].decode("ascii")
    tcout = [l.split() for l in tcout.split("\n")[2:-1]] #assumes first two lines and last line are not data (FIXME: fragile code)
    print("-- Trace lines are : "+", ".join([data(strip(l[0])) for l in tcout]))
    mineff = min([int(l[1]) for l in tcout])
    maxrel = max([float(l[2]) for l in tcout])
    print("-- Minimum effsize is "+data(mineff)+" and maximum rel_diff is "+data(maxrel))
    if mineff < 50 or maxrel > 0.3:
        print(bad("The chains do not seem to have converged!"))
    else:
        print(good("The chains seem to have converged!"))
        convergence = True
else:
    print("-- "+boldred("Warning")+" Something went wrong with tracecomp!\n-- Output of command was:\n")
    print(p4.communicate()[1].decode("ascii"))

#===================================================================================================
print(step("Convergent site analysis"))
detection = False

print("-- Running readdiffsel")

def rdrun(no):
    p5 = Popen(["_build/readdiffsel -x "+str(burnin)+" 1 "+str(iterations)+" "+chainname+no], shell=True, stdout=PIPE, stderr=PIPE)
    r5 = p5.wait()

    meanname = chainname+no+"_1.meandiffsel"
    signname = chainname+no+"_1.signdiffsel"

    if r5 == 0 and isfile(meanname) and isfile(signname):
        print("-- Readdiffsel executed; results stored in "+data(meanname)+" and "+data(signname))
    else:
        print(boldred("ERROR")) # TODO
        print(p5.communicate()[1].decode("ascii"))
        exit(1)

    rdoutfile = open(signname, 'r')
    rdout = [l.split() for l in rdoutfile.readlines()]
    aas = [int(l[1]) for l in rdout]
    print("-- Convergence detected for aminoacids "+", ".join([data(aa) for aa in aas]))
    return aas

if rdrun("") == [11] and rdrun("2") == [11]:
    print(good("The chain detected the correct aminoacids!"))
    detection = True
else:
    print(bad("The chains detected the wrong aminoacids!"))

#===================================================================================================
print(step("Non-regression test result"))
if detection and convergence:
    print(success("Test has succeeded!"))
else:
    print(failure("Test failed!"))
    if not allowfailure:
        exit(1)
