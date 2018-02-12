# Copyright or Copr. Centre National de la Recherche Scientifique (CNRS) (2017/11/27)
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


from ete3 import*
from diffsel_script_utils import *

#===================================================================================================
print(step("Parsing command line arguments"))

from argparse import ArgumentParser, FileType
parser = ArgumentParser(description='Annotates a phylogenetic trees with condition numbers compatible with diffsel.')
parser.add_argument('inputFile', metavar="input", type=FileType('r'), nargs=1, help='the sequence file (phylip format)')
parser.add_argument('-s', '--sister-branch-cond', dest="sister", action='store_true', help="toggle the use of a different condition for the sister branches of convergent branches")

args = parser.parse_args()

tree_file = args.inputFile[0]
print("-- Sequence file is "+param(tree_file.name))
sister = args.sister
print("-- Sister branch condition: "+param(sister))
out_file = tree_file.name+".annotated"
print("-- Output file is "+param(out_file))


#===================================================================================================
print(step("Tree retrieval and preparation"))

print("-- Reading tree from file")
t = Tree(tree_file.name)

print("-- Setting all branch lengths to "+data(0))
for n in t.traverse():
    n.dist = 0 # using dist as condition


#===================================================================================================
print(step("Convergent branch selection"))

nstyle = NodeStyle()
nstyle["size"] = 20
nstyle["fgcolor"] = "darkred"

nstyle2 = NodeStyle()
nstyle2["size"] = 20
nstyle2["fgcolor"] = "blue"

nstyle3 = NodeStyle()
nstyle3["size"] = 20
nstyle3["fgcolor"] = "green"

print("-- Numbering nodes")
i=0
for n in t.traverse():
    label = TextFace(str(i))
    n.add_face(label, column=0, position="branch-top")
    n.set_style(nstyle2)
    i += 1


print("-- Starting subtree selection")
while True:
    t.show()
    data = raw_input(ask_input("Please enter start of convergent subtree: "))
    if data.isdigit():
        nb = int(data)
        i = 0
        for st_canditate in t.traverse():
            if i == nb:
                for n in st_canditate.traverse():
                    n.dist=1
                    n.set_style(nstyle)
                if sister:
                    for n in st_canditate.get_sisters()[0].traverse():
                        n.dist=2
                        n.set_style(nstyle3)
            i += 1
    else:
        print("-- Input was not an integer; finished selection of subtrees")
        break

#===================================================================================================
print(step("Writing result to file"))
t.write(format=1, outfile = out_file)
