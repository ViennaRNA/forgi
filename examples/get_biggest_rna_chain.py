#!/usr/bin/python

import sys
from optparse import OptionParser

import forgi.threedee.utilities.pdb as ftup

def main():
    usage = """
    python get_biggest_rna_chain.py in.pdb out.pdb

    Extract the largest RNA chain from the file in.pdb and store it in out.pdb.
    """
    num_args= 2
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    c = ftup.get_biggest_chain(args[0])
    ftup.output_chain(c, args[1])

if __name__ == '__main__':
    main()

