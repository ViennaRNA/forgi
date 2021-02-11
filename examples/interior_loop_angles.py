#!/usr/bin/python

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud

import os.path as op
import sys
from optparse import OptionParser
import csv


def main():
    usage = """
    python interior_loop_angles.py pdb_file

    Iterate over the interior loop angles and calculate how much of a kink
    they introduce between the two adjacent stems.
    """
    num_args = 0
    parser = OptionParser(usage=usage)

    parser.add_option("-o","--output", action="store", help="Store data in csv with this filename")
    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    data = list()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    cg, = ftmc.CoarseGrainRNA.from_pdb(op.expanduser(args[0]))
    for iloop in cg.iloop_iterator():
        conn = cg.connections(iloop)
        angle = ftuv.vec_angle(
            cg.coords[conn[0]][1] - cg.coords[conn[0]][0], cg.coords[conn[1]][1] - cg.coords[conn[1]][0])
        data.append([iloop,angle])

        fud.pv('iloop, angle')

    if options.output and len(data):
        with open(options.output, 'w') as FILE:
            writer=csv.writer(FILE,delimiter="\t", lineterminator="\n")
            writer.writerow(["iloop", "angle"])
            for row in data:
                writer.writerow(row)

if __name__ == '__main__':
    main()
