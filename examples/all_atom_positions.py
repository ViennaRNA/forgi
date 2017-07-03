#!/usr/bin/python

from __future__ import print_function
from builtins import map
import collections as c
import itertools as it

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.vector as ftuv

import forgi.utilities.debug as fud
import numpy as np
import sys
from optparse import OptionParser

def main():
    usage = """
    python average_atom_positions.py file1.pdb file2.pdb ...
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')
    parser.add_option('-p', '--pseudoknots', dest='pseudoknots', default=False, action='store_true', help='Allow pseudoknots in the CG structure')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    poss = c.defaultdict(list)
    sources = c.defaultdict(list)

    for i,arg in enumerate(args):
        cg = ftmc.from_pdb(arg, remove_pseudoknots=not options.pseudoknots)

        if len(list(cg.stem_iterator())) == 0:
            print("skipping {}: no stems".format(arg), file=sys.stderr)
            continue

        for d in cg.defines.keys():
            if np.allclose(cg.coords[d][0], cg.coords[d][1]):
                print("File {}: Degenerate coordinates for element: {}".format(i, d), file=sys.stderr)
                continue

            origin, basis = ftug.element_coord_system(cg, d)

            if d[0] == 'i' or d[0] == 'm':
                conn = cg.connections(d)
                conn_type = cg.connection_type(d, conn)
            else:
                conn_type = 0

            for i, r in it.izip(it.count(),
                                cg.define_residue_num_iterator(d)):

                # add only the base atoms which are relevant to the calculation
                # of the chi torsion angle
                resname = cg.chain[cg.seq_ids[r-1]].resname.strip()

                if resname not in ftup.chi_torsion_atoms.keys():
                    print("Unknown nucleotide name:", resname, file=sys.stderr)
                    continue

                atoms = ftup.nonsidechain_atoms + ftup.chi_torsion_atoms[resname][-2:]
                scatoms=ftup.side_chain_atoms[resname]
                for aname in atoms+scatoms:
                    try:
                        a = cg.chain[cg.seq_ids[r-1]][aname]
                    except KeyError as ke:
                        # missing an atom
                        continue

                    # The C1'->B1 and B1->B2 vectors define the plane of the base
                    # The O4'->C1'->B1->B2 sequence defines the torsion 
                    # angle chi
                    if aname == ftup.chi_torsion_atoms[resname][-2]:
                        aname = 'B1'
                    elif aname == ftup.chi_torsion_atoms[resname][-1]:
                        aname = 'B2'
                    elif aname in scatoms:
                        aname=resname+"."+aname
                    avec = a.get_vector().get_array()
                    atom_pos = ftuv.change_basis(avec - origin, basis, ftuv.standard_basis)
                    identifier = "%s %s %d %d %s" % (d[0], 
                                                  " ".join(map(str, cg.get_node_dimensions(d))),
                                                  conn_type,                              
                                                  i, aname)
                    poss[identifier] += [atom_pos]
                    sources[identifier] += [d]

                    print("{}:{}".format(identifier, ",".join(map(str, atom_pos))))

if __name__ == '__main__':
    main()

