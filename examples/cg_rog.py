#!/usr/bin/python

import sys
from optparse import OptionParser

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.model.descriptors as ftmd

def main():
    usage = """
    python cg_rmsd.py file1.cg

    Calculate the ROG of a coarse grain models.
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    cg1 = ftmc.CoarseGrainRNA(args[0])


    coords = cg1.get_ordered_stem_poss()
    print ftmd.radius_of_gyration(coords)
if __name__ == '__main__':
    main()

