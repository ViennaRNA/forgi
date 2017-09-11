#!/usr/bin/python

from __future__ import print_function
from builtins import zip
from builtins import range
import forgi.graph.bulge_graph as fgb
import forgi.utilities.debug as fud

import sys
from optparse import OptionParser

def main():
    usage = """
    python fasta_to_reduced_graph.py fasta_file

    Generate a reduced graph representation of this secondary structure.
    fasta_file should be either a filename pointing to a fasta file
    or a dash (-) indicating that the input should come from stdin.
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    if args[0] == '-':
        text = sys.stdin.read()
    else:
        with open(args[0], 'r') as f:
            text = f.read()

    # load a BulgeGraph from a fasta-like file
    # i.e.
    # >test
    # AACCGG
    # ((..))
    bg = fgb.from_fasta_text(text)

    prev=None
    seen = set()
    out_str = ''
    elements_traversed = []

    for i in range(1, bg.seq_length + 1):
        # iterate over each nucleotide and get the name of the element
        # that it is part of
        node = bg.get_node_from_residue_num(i)

        if node == prev:
            # we've already seen this element
            continue
        if node[0] == 's':
            if node in seen:
                # we've seen this stem before, so we just need to close the bracket
                out_str += ')'
            else:
                # new stem, so open a bracket
                out_str += '('

            elements_traversed += [node]
        elif node[0] == 'i':
            out_str += '.'
            elements_traversed += [node]

        prev = node
        seen.add(node)

    print(out_str)
    print(",".join(elements_traversed))

    for t,e in zip(out_str, elements_traversed):
        # print the sequences of each element in the reduced representation
        print(t, e, bg.get_define_seq_str(e))

if __name__ == '__main__':
    main()

