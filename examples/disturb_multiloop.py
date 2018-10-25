#!/usr/bin/python

from __future__ import print_function
import sys
import random
import os.path as op
import forgi.graph.bulge_graph as cgb
import forgi.threedee.model.coarse_grain as cmc
import forgi.utilities.commandline_utils as fuc
from optparse import OptionParser
import logging


def main(args):

    with fuc.hide_traceback():
        bg, = fuc.cgs_from_args(args, 1, "any", enable_logging=True)

    multiloops, _ = bg.find_multiloop_loops()
    for multi in multiloops:
        shortened = set()
        for m in multi:
            if m[0] != 'm':
                continue
            connected_stems = list(bg.edges[m])
            for to_shorten in connected_stems:
                if to_shorten in shortened:
                    continue
                shortened.add(to_shorten)

                # find the stems which are connected to this multiloop
                # and pick a random one
                db = list(bg.to_dotbracket_string())

                # get the side of the stem which is connected to the
                # multiloop
                (s1b, s1e) = bg.get_sides(to_shorten, m)

                # print to_shorten, s1b, "(", bg.defines[to_shorten], ")"
                # the nucleotides that need to be changed
                to_change = bg.get_side_nucleotides(to_shorten, s1b)
                # print bg.defines[to_shorten], to_change

                db[to_change[0] - 1] = '.'
                db[to_change[1] - 1] = '.'
                print("".join(db))


parser = fuc.get_rna_input_parser("Add (default) or remove a basepair from a random multiloop in the file and "
                                  "print out the resulting structure.", nargs=1, rna_type="any", enable_logging=True)
args = parser.parse_args()
if __name__ == '__main__':
    main(args)
