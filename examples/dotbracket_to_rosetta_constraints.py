#!/usr/bin/python

import sys, math
import forgi.graph.bulge_graph as cgb

from optparse import OptionParser

def print_rosetta_constraints(bg):
    '''
    Convert the stems into base pair constraints compatible
    with the rna_denovo program.

    @param bg: A BulgeGraph data structure
    '''
    for s in bg.stem_iterator():
        for i in range(bg.stem_length(s)):
            print "STEM PAIR %d %d" % (bg.defines[s][0] + i, bg.defines[s][3] - i)

def main():
    usage = """
        Usage: ./dotbracket_to_rosetta_constraints.py

        Convert a dotbracket file to a set of constraints for the Rosetta rna_denovo program.
        These constraints will describe the stems within the structure by 
        showing which necleotides need to pair with which other nucleotides.
        """
    parser = OptionParser(usage = usage)
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)
    if args[0] == '-':
        f = sys.stdin
    else:
        f = open(args[0])

    brackets = "".join(f.readlines()).replace('\n', '')
    bg = cgb.BulgeGraph()
    bg.from_dotbracket(brackets)
    print_rosetta_constraints(bg)

if __name__=="__main__":
    main()
