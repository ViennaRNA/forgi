#!/usr/bin/python

import sys
import argparse

import forgi.threedee.model.coarse_grain as ftmc
import forgi.projection.projection2d as ftmp
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.rmsd as ftur

import numpy as np



def get_parser():
    parser = argparse.ArgumentParser(description='Calculate the RMSD in 2D space between 2 projections.')
    parser.add_argument('files', type=str, nargs=2, help="The *.coord files")
    parser.add_argument('--directions', nargs=2, type=str, help="The projection directions: e.g. 1.0,1.0,2.4 2.4,0,1")
    parser.add_argument('--plot', type=bool, default=False help="Plot the projections")
    return parser


def main(parser):
    args=parser.parse_args()

    cg1 = ftmc.CoarseGrainRNA(args.files[0])
    cg2 = ftmc.CoarseGrainRNA(args.files[1])
  
    dir1=np.array(args.directions[0].split(","), dtype=float)
    dir2=np.array(args.directions[1].split(","), dtype=float)
    
    proj1=ftmp.Projection2D(cg1, dir1)
    proj2=ftmp.Projection2D(cg2, dir2)

    vrs1 = np.array([x for p in sorted(proj1._coords.keys()) for x in proj1._coords[p]])
    vrs2 = np.array([x for p in sorted(proj2._coords.keys()) for x in proj2._coords[p]])

    print ftur.centered_rmsd(vrs1, vrs2)
    if args.plot:
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots()
        proj1.plot(ax, line2dproperties={"color":"green"})
        proj2.plot(ax, line2dproperties={"color":"red"})
        plt.show()


parser=get_parser()
if __name__ == '__main__':
    main(parser)

