#!/usr/bin/python

from __future__ import print_function
import forgi.graph.bulge_graph as fgb
import sys
from optparse import OptionParser

def main():
    usage = """
    python cg_to_ss_fasta.py file.cg

    Convert a coarse-grain RNA file to a dotbracket string.
    BEWARE: Structures with pseudoknots will give bogus output.
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    cg = fgb.BulgeGraph(args[0])
    print(cg.to_fasta_string())

if __name__ == '__main__':
    main()

