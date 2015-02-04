#!/usr/bin/python

import Bio.PDB as bpdb
import Bio.KDTree as bkd

import forgi.threedee.utilities.pdb as ftup
import numpy as np

import sys
from optparse import OptionParser

def main():
    usage = """
    python burial.py pdb_file

    Calculate how buried each nucleotide is in this PDB file.
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    chain = ftup.load_structure(args[0])
    atoms = bpdb.Selection.unfold_entities(chain, 'A')
    atom_coords = np.array([a.get_coord() for a in atoms])

    kd = bkd.KDTree(3)
    kd.set_coords(atom_coords)

    for r in chain.get_list():
        pass


if __name__ == '__main__':
    main()

