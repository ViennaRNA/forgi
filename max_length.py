import Bio.PDB as bpdb
import sys
import itertools as it


import forgi.graph.bulge_graph as cgb
import forgi.threedee.model.coarse_grain as cmc

# load the pdb structure
s = bpdb.PDBParser().get_structure('blah', sys.argv[1])

# load a bulge-graph representation
m = cmc.from_pdb(sys.argv[1])
bg = m.bg


# little function to see if a given residue number corresponds
# to a particular named element of the bulge graph
def is_inelement(bg, element, resnum):
    d = bg.defines[element]

    for i in range(0, len(d), 2):
        if resnum > d[i] and resnum < d[i+1]:
            return True

    return False

# take only the C1' atom of the structure and make sure it's not in
# the fiveprime or threeprime element
atoms = [a for a in bpdb.Selection.unfold_entities(s, 'A') if (a.name == "C1'" and not is_inelement(bg, 'f1', a.parent.id[1]) and not is_inelement(bg, 't1', a.parent.id[1]))]

# calculate all the distances and return the maximum
dists = [(a1 - a2, a1.parent, a2.parent) for a1, a2 in it.combinations(atoms, 2)]
#print dists
print max(dists)

