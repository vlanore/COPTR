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

# thanks to http://www.petercollingridge.co.uk/book/export/html/474
bases = ['a', 'c', 'g', 't']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

reverse_bases = dict(zip(bases, range(4)))
short_aa_table = "FFLLSSSSYYCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
# short_aa_table = "FFLLSSSSYYCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
aas = "".join(sorted("FSYCWLPHQRIMTNKVADEG"))
reverse_aas = dict(zip(aas, range(20)))

def cno_to_nuc(cno):
    return codons[cno]

# def cno_to_nuc(cno):
#     if (cno <10):
#         return codons[cno]
#     elif (cno == 10 or cno == 11):
#         return codons[cno+1]
#     else:
#         return codons[cno+2]

def synonym(c1, c2):
    return codon_table[cno_to_nuc(c1)] == codon_table[cno_to_nuc(c2)]

def diff(cno1, cno2):
    c1 = cno_to_nuc(cno1)
    c2 = cno_to_nuc(cno2)
    r = []
    for i in range(len(c1)):
        if c1[i] != c2[i]:
            r.append((c1[i], c2[i]))
    return r

def case(c1, c2):
    d = diff(c1, c2)
    if (len(d) != 1):
        pass
        # print("\tQcodons[cod][%d][%d] := 0" % (c1, c2))
    elif synonym(c1, c2):
        print("\tQcodons[cod][%d][%d] := Q[%d][%d]" % (c1+1, c2+1, reverse_bases[d[0][0]]+1, reverse_bases[d[0][1]]+1))
    elif (amino_acids[c1] == '*' or amino_acids[c2] == '*'):
        pass
        # print("\tQcodons[cod][%d][%d] := 0" % (c1, c2))
    else:
        print("\tQcodons[cod][%d][%d] := Q[%d][%d] * sqrt(fitness[cod][%d]/fitness[cod][%d])" %
              (c1+1, c2+1, reverse_bases[d[0][0]]+1, reverse_bases[d[0][1]]+1, reverse_aas[amino_acids[c2]]+1, reverse_aas[amino_acids[c1]]+1))


def print_all():
    print("for (cod in 1:nsites) {")
    print("  for (i in 1:64) {")
    print("    for (j in 1:64) {")
    print("      Qcodons[cod][i][j] := abs(0)")
    print("    }")
    print("  }")
    print("}")
    print("for (cod in 1:nsites) {")
    for cno in range(64):
        for cno2 in range(64):
            case(cno, cno2)
    print("\tR[cod] := fnFreeK(Qcodons[cod])")
    print("}")

print_all()
