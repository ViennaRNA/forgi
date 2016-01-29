#!/usr/bin/python

import sys
import argparse

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.projection2d as ftmp
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.rmsd as ftur

import numpy as np
import pickle

from collections import OrderedDict

def get_parser():
    parser = argparse.ArgumentParser(description='For different offsets, calculates the average  RMSD in the optimal '
                         'projection and in 3D space for all pairs of structures with given offsets.')
    parser.add_argument('files', type=str, nargs='*', help="An ORDERED list of *.coord files")
    parser.add_argument('--step-size', type=int, help="How many steps two subsequent files are apart. Used for the x-axis of the plot.", default=1)
    parser.add_argument('--max-diff', type=str, help="Maximal offset in FILES to be considered, or a ','-separated list of offsets to be considered.", default=20)
    parser.add_argument('--target-structure', type=str, help="In addition, calculate the RMSD of every step to structure given here.")
    parser.add_argument('--save-fig', type=str, help="Save the figure (by pickleing the 'matplotlib.figure.Figure' object) to the given filename. "
                                                     "WARNING: Pickle is unsave. Dont load the file if someone else created it!")
    return parser


def main(args):

    cgs=[]
    projs=[]
    for file_ in args.files:
        cgs.append(ftmc.CoarseGrainRNA(file_))
        projs.append(ftmp.Projection2D(cgs[-1]))

    p_rmsds=OrderedDict()
    cg_rmsds=OrderedDict()
    p_rmsds[0]=0.
    cg_rmsds[0]=0.

    if "," in args.max_diff:
        diffs=map(int,args.max_diff.split(","))
    else:
        diffs=range(1,int(args.max_diff)+1)
    for diff in diffs:
        print("diff {}".format(diff))
        prmsd=0
        cgrmsd=0
        count=0
        for i in range(0,len(cgs)-diff):
            count+=1
            vrs1 = np.array([x for p in sorted(projs[i]._coords.keys()) for x in projs[i]._coords[p]])
            vrs2 = np.array([x for p in sorted(projs[i+diff]._coords.keys()) for x in projs[i+diff]._coords[p]])
            prmsd+=ftur.centered_rmsd(vrs1, vrs2)
            vrs1 = ftug.bg_virtual_residues(cgs[i])
            vrs2 = ftug.bg_virtual_residues(cgs[i+diff])
            cgrmsd+=ftur.centered_rmsd(vrs1,vrs2)
        if count:
            p_rmsds[diff]=prmsd/count
            cg_rmsds[diff]=cgrmsd/count
    print "projection RMSDs:", p_rmsds
    print "3D-RMSDs:", cg_rmsds
    import matplotlib.pyplot as plt
    if args.target_structure:
        target=ftmc.CoarseGrainRNA(args.target_structure)
        target_rmsds=[]
        target_proj_rmsds=[]
        try:
            target_proj=ftmp.Projection2D(target)
        except ValueError:
            pass
        else:                
            target_vrs = np.array([x for p in sorted(target_proj._coords.keys()) for x in target_proj._coords[p]])
            for proj in projs:
                vrs1 = np.array([x for p in sorted(proj._coords.keys()) for x in proj._coords[p]])
                prmsd=ftur.centered_rmsd(vrs1, target_vrs)
                target_proj_rmsds.append(prmsd)            
        target_vrs = ftug.bg_virtual_residues(target)
        for cg in cgs:
            vrs1 = ftug.bg_virtual_residues(cg)
            cgrmsd=ftur.centered_rmsd(vrs1,target_vrs)
            target_rmsds.append(cgrmsd)
        fig,(ax,ax2)=plt.subplots(2)
        xval=np.arange(len(cgs))*args.step_size
        if target_proj_rmsds:
            ax2.plot(xval,target_proj_rmsds,label="Projection", color="green")
        ax2.plot(xval,target_rmsds,label="3D", color="blue")
        ax2.set_title("RMSD to target structure")
        ax2.set_xlabel("step")
        ax2.set_ylabel("RMSD")
        leg=ax2.legend(framealpha=0.8, fancybox=True)
        leg.get_frame().set_linewidth(0.0)
        plt.subplots_adjust(hspace=0.4)
    else:
        fig,ax=plt.subplots()
    p_rmsds=np.array(p_rmsds.items())
    cg_rmsds=np.array(cg_rmsds.items())
    ax.plot(p_rmsds[:,0]*args.step_size, p_rmsds[:,1], label="Projection",color="green")
    ax.plot(cg_rmsds[:,0]*args.step_size, cg_rmsds[:,1], label="3D", color="blue")
    ax.set_title("Average RMSD between structures X steps apart.")
    ax.set_xlabel("steps apart")
    ax.set_ylabel("RMSD")
    leg=ax.legend(framealpha=0.8, fancybox=True)
    leg.get_frame().set_linewidth(0.0)
    fig.patch.set_facecolor('white')
    if args.save_fig:
        with open(args.save_fig, "w") as f:
            pickle.dump(fig, f)
    plt.show()

parser=get_parser()
if __name__ == '__main__':    
    args=parser.parse_args()
    main(args)
