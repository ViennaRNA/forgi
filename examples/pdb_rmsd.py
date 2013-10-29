#!/usr/bin/python

import sys
import corgy.utilities.pdb as cup
from optparse import OptionParser

def main():
    usage = """
    ./pdb_rmsd.py pdb_file1 pdb_file2

    Extract the largest RNA chain from each file and calculate the rmsd between
    all the matching atoms. Matching atoms are those that share a position and name.

    I.e. residue 1, atom C1'
    """
    num_args= 2
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    chain1 = cup.load_structure(args[0])
    chain2 = cup.load_structure(args[1])

    print cup.pdb_rmsd(chain1, chain2)

if __name__ == '__main__':
    main()

