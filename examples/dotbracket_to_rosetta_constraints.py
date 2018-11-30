#!/usr/bin/python

from __future__ import print_function
from builtins import range
import sys
import math
import forgi.graph.bulge_graph as cgb
import forgi.utilities.commanline_utils as fuc
from optparse import OptionParser


def print_rosetta_constraints(bg):
    '''
    Convert the stems into base pair constraints compatible
    with the rna_denovo program.

    @param bg: A BulgeGraph data structure
    '''
    for s in bg.stem_iterator():
        for i in range(bg.stem_length(s)):
            print("STEM PAIR %d %d" %
                  (bg.defines[s][0] + i, bg.defines[s][3] - i))


def main():
    usage = """
        Convert a dotbracket file to a set of constraints for the Rosetta rna_denovo program.
        These constraints will describe the stems within the structure by
        showing which necleotides need to pair with which other nucleotides.
        """
    parser = fuc.get_rna_input_parser(usage, nargs=1)
    args = parser.parse_args()

    cgs = fuc.cgs_from_args(args)
    if len(cgs)>1:
        raise ValueError("More than one RNA molecule found in the input.")
    print_rosetta_constraints(cgs[0])


if __name__ == "__main__":
    main()
