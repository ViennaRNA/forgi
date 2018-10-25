#!/usr/bin/python

from __future__ import print_function
import sys
import argparse

import forgi.utilities.commandline_utils as fuc
import forgi.threedee.model.coarse_grain as ftmc
import forgi.projection.projection2d as ftmp
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.model.similarity as ftms

import numpy as np


def get_parser():
    parser = fuc.get_rna_input_parser('Calculate the RMSD in 2D space between 2 projections of '
                                      'coarse-grained RNAs.', nargs=2, rna_type="3d",
                                      enable_logging=True)
    parser.add_argument('--directions', nargs=2, type=str,
                        help="The projection directions: e.g. 1.0,1.0,2.4 2.4,0,1", required=True)
    parser.add_argument('--plot', action="store_true",
                        help="Open a window with a plot of the projections")
    return parser


def main(parser):
    args = parser.parse_args()

    with fuc.hide_traceback():
        cg1, cg2 = fuc.cgs_from_args(
            args, nargs=2, rna_type="3d", enable_logging=True)

    dir1 = np.array(args.directions[0].split(","), dtype=float)
    dir2 = np.array(args.directions[1].split(","), dtype=float)

    proj1 = ftmp.Projection2D(cg1, dir1)
    proj2 = ftmp.Projection2D(cg2, dir2)

    vrs1 = np.array([x for p in sorted(proj1._coords.keys())
                     for x in proj1._coords[p]])
    vrs2 = np.array([x for p in sorted(proj2._coords.keys())
                     for x in proj2._coords[p]])

    print(ftms.rmsd(vrs1, vrs2))
    if args.plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        proj1.plot(ax, line2dproperties={"color": "green"})
        proj2.plot(ax, line2dproperties={"color": "red"})
        plt.show()


parser = get_parser()
if __name__ == '__main__':
    main(parser)
