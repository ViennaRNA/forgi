#!/usr/bin/python

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.vector as ftuv
import itertools as it

import sys
from optparse import OptionParser

def main():
    usage = """
    usage
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    cg = ftmc.from_pdb(args[0])

    angles = []
    for loop in it.chain(cg.iloop_iterator(), cg.mloop_iterator()):
        conn = cg.connections(loop)

        (s1b, s1e) = cg.get_sides(conn[0], loop)
        (s2b, s2e) = cg.get_sides(conn[1], loop)

        angle = ftuv.vec_angle(cg.coords[conn[0]][s1b] - cg.coords[conn[0]][s1e], 
                               cg.coords[conn[1]][s2e] - cg.coords[conn[1]][s2b])

        for rn in cg.define_residue_num_iterator(loop, adjacent=True):
            angles += [(rn, angle)]

    for rn, angle in sorted(angles):
        print "{}:{}".format(rn, angle)

if __name__ == '__main__':
    main()

