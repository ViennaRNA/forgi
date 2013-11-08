#!/usr/bin/python

import sys
import forgi.model.coarse_grain as cmc
import os.path as op
from optparse import OptionParser

def main():
    usage = """
    ./bulge_graph_from_pdb.py pdb_file 

    Create a bulge graph from the information contained in the tertiary
    structure of an RNA molecule.
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)
    
    cg = cmc.from_pdb(args[0])

if __name__ == '__main__':
    main()

