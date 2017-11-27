#!/usr/bin/python3.5

# Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2017/05/03)
# Contributors:
# - Vincent Lanore <vincent.lanore@gmail.com>

# This software is a computer program whose purpose is to provide the necessary classes to write ligntweight component-based
# c++ applications.

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

#===================================================================================================
print(step("Parsing command line arguments"))

from argparse import ArgumentParser, FileType
parser = ArgumentParser(description='Removes indels from a sequence file.')
parser.add_argument('inputFile', metavar="input", type=FileType('r'), nargs=1, help='the sequence file (phylip format)')

args = parser.parse_args()

seq_file = args.inputFile[0]
print("-- Sequence file is "+param(seq_file.name))
out_file = seq_file.name+".woindels"
print("-- Output file is "+param(out_file))

#===================================================================================================
print(step("Extracting and verifying data from file"))

nb_taxa, length = [int(x) for x in seq_file.readline().split()]
print("-- Number of taxa: "+data(nb_taxa))
print("-- Sequence length before indel removal: "+data(length))

sequences = [l.split() for l in seq_file]
if (nb_taxa == len(sequences)):
    print(good("Number of sequences matches number of taxa"))
else:
    print(failure("number of sequences is "+str(len(sequences))+" instead of "+str(nb_taxa)))
    exit(1)
for i in range(nb_taxa):
    if (len(sequences[i][1]) != length):
        print(failure("sequence "+str(i)+" ("+sequences[i][0]+") has length "+str(len(sequences[i][1]))+" instead of "+str(length)))
        exit(1)
print(good("Lengths of sequences all match expected sequence length"))

#===================================================================================================
print(step("Removing indels!"))
indel_pos = []
for s in sequences:
    for pos in range(length):
        if s[1][pos] == '-' or s[1][pos] == '?':
            if not (pos in indel_pos):
                indel_pos.append(pos)
print("-- Found "+data(len(indel_pos))+" positions with indels")
print("-- Building resulting alignment... ", end="")
non_indel_pos = set(range(length)).difference(indel_pos)
result = [[s[0], [s[1][i] for i in non_indel_pos]] for s in sequences]
print("Done")
print("-- Resulting alignment length is "+data(len(result[1][1]))+" ("+data(round(100.*len(indel_pos)/length))+"% of sites were removed)")

#===================================================================================================
print(step("Writing result to file"))
out = open(out_file, 'w')
out.write(str(nb_taxa)+'\t'+str(len(non_indel_pos))+'\n')
for s in result:
    out.write(s[0]+'\t')
    out.write(str(''.join(s[1])))
    out.write('\n')
out.close()
