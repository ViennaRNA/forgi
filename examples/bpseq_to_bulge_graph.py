#!/usr/bin/python

import forgi.graph.bulge_graph as fgb

import sys
from optparse import OptionParser

def main():
    usage = """
    python bpseq_to_bulge_graph.py secondary_structure.bpseq
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    with open(args[0], 'r') as f:
        text = f.read()
        bg = fgb.BulgeGraph()
        bg.from_bpseq_str(text)
        print bg.to_bg_string()

if __name__ == '__main__':
    main()

