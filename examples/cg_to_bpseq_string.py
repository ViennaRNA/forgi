#!/usr/bin/python

import forgi.graph.bulge_graph as fgb
import sys
from optparse import OptionParser

def main():
    usage = """
    python cg_to_bp_string.py file.cg

    Convert a coarse-grain file to a bpseq string.
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
    print cg.to_bpseq_string()

if __name__ == '__main__':
    main()

