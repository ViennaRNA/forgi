#!/usr/bin/python

import sys
from optparse import OptionParser

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.rmsd as ftur

def main():
    usage = """
    python cg_rmsd.py file1.cg file2.cg

    Calculate the RMSD between two coarse grain models.
    """
    num_args= 2
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    cg1 = ftmc.CoarseGrainRNA(args[0])
    cg2 = ftmc.CoarseGrainRNA(args[1])

    vrs1 = ftug.bg_virtual_residues(cg1)
    vrs2 = ftug.bg_virtual_residues(cg2)

    print ftur.centered_rmsd(vrs1, vrs2)

if __name__ == '__main__':
    main()

