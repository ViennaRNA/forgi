#!/usr/bin/python

import sys
import corgy.graph.bulge_graph as cgb

from optparse import OptionParser

def main():
    usage = """
    ./longest_stem.py dotbracket_file
    """
    num_args=1
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    if args[0] == '-':
        f = sys.stdin
    else:
        f = open(args[0])

    brackets = "".join(f.readlines()).replace('\n', '')
    bg = cgb.BulgeGraph()
    bg.from_dotbracket(brackets)

    biggest_stem = (-1, 'x')

    for s in bg.stem_iterator():
        if bg.stem_length(s) > biggest_stem[0]:
            biggest_stem = (bg.stem_length(s), s)

    print biggest_stem[0]

if __name__ == '__main__':
    main()

