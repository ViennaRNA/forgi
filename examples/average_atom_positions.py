#!/usr/bin/python

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
    parser.add_option('-a', '--all', dest='all_entries', default=False, action='store_true', help='Store all positions')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    poss = c.defaultdict(list)

    for i,arg in enumerate(args):
        cg = ftmc.from_pdb(arg)
        fud.pv('i')

        for d in cg.defines.keys():
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
                resname = cg.chain[r].resname.strip()
                atoms = ftup.nonsidechain_atoms + ftup.chi_torsion_atoms[resname][-2:]

                for aname in atoms:
                    try:
                        a = cg.chain[r][aname]
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

                    avec = a.get_vector().get_array()
                    atom_pos = ftuv.change_basis(avec - origin, basis, ftuv.standard_basis)
                    identifier = "%s %s %d %d %s" % (d[0], 
                                                  " ".join(map(str, cg.get_node_dimensions(d))),
                                                  conn_type,                              
                                                  i, aname)
                    poss[identifier] += [atom_pos]

    print "import collections as co"

    if options.all_entries:
        print "all_atom_poss = dict()"
        for key in poss.keys():
            print 'all_atom_poss["%s"] = [%s] #%d' % (key, 
                                                      ",".join(["[%s]" % (",".join(map(str, pos))) for pos in poss[key]]), len(poss[key]))
    else:
        print "avg_atom_poss = dict()"
        for key in poss.keys():
            pos = np.mean(poss[key], axis=0)
            print 'avg_atom_poss["%s"] = [%s] #%d' % (key, ",".join(map(str, pos)), len(poss[key]))

if __name__ == '__main__':
    main()

