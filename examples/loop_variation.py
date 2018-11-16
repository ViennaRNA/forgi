#!/usr/bin/python

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.stats as ftms
import forgi.utilities.debug as fud

import itertools as it

import os.path as op
import sys
from optparse import OptionParser


def calculate_variation(angle_stats, loop_size):
    '''
    Calculate how much the statistics for a particular vary
    based on its dimensions. I.e., small bulges should vary
    very little whereas larger loops should vary a lot.

    An exact method for calculating how much a loop can vary
    is difficult to derive, but this method will use the
    simplest available, the volume of the n-dimensional
    enclosure defined by the minimum and the maximum coordinates.

    :param stats: forgi.threedee.model.stats.AngleStats
    :param dims: The dimensions of the loop (i.e. (1,3))
    :return: The volume of the accessible area.
    '''
    # ang_type indicates whether it's an iloop forward/backward
    # or a multiloop forward/backward
    ang_types = [1, 2, 3, 4]

    for ang_type in ang_types:
        ang_dims = tuple(list(loop_size) + [ang_type])
        if ang_dims in angle_stats:

            fud.pv('ang_dims')
            fud.pv('len(angle_stats[ang_dims])')


def main():
    usage = """
    python loop_variation.py pdb_file|cg_file

    Calculate how much the geometry of each loop-associated nucleotide can
    vary based on the size of the loop that it is in.

    Whether the input is a pdb or a cg file depends on the extension.
    """
    num_args = 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    fn, fe = op.splitext(args[0])

    if fe == '.cg':
        cg = ftmc.CoarseGrainRNA(args[0])
    else:
        cg = ftmc.from_pdb(args[0])

    angle_stats = ftms.get_angle_stats()

    for loop in it.chain(cg.iloop_iterator(), cg.mloop_iterator()):
        calculate_variation(angle_stats, cg.get_bulge_dimensions(loop))
        #fud.pv('loop, cg.get_bulge_dimensions(loop), cg.get_angle_type(loop)')


if __name__ == '__main__':
    main()
