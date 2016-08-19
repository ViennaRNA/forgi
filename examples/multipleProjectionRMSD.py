#!/usr/bin/python

import sys
import argparse

import forgi.threedee.model.coarse_grain as ftmc
import forgi.projection.projection2d as ftmp
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.rmsd as ftur
import forgi.threedee.model.similarity as ftuc
import forgi.projection.hausdorff as fph

import numpy as np
import pickle, math
import scipy.ndimage
import itertools as it

from collections import OrderedDict

def get_parser():
    parser = argparse.ArgumentParser(description='For different offsets, calculates the average  RMSD in the optimal '
                         'projection and in 3D space for all pairs of structures with given offsets.')
    parser.add_argument('files', type=str, nargs='*', help="An ORDERED list of *.coord files")
    parser.add_argument('--step-size', type=int, help="How many steps two subsequent files are apart. Used for the x-axis of the plot.", default=1)
    parser.add_argument('--max-diff', type=str, help="Maximal offset in FILES to be considered, or a ','-separated list of offsets to be considered.", default="20")
    parser.add_argument('--target-structure', type=str, help="In addition, calculate the RMSD of every step to structure given here.")
    parser.add_argument('--hausdorff-reference', type=str, help="In addition, calculate the Hausdorffdistance of every step to the png image file given here.")
    parser.add_argument('--hausdorff-scale', type=float, help="Side length of the hausdorff-reference image in Angstrom.")
    parser.add_argument('--save-fig', type=str, help="Save the figure (by pickleing the 'matplotlib.figure.Figure' object) to the given filename. "
                                                     "WARNING: Pickle is unsave. Dont load the file if someone else created it!")
    return parser

def get_refimg_longest_axis(ref_img):
    max_sdist=0
    best_points=None
    for p1, p2 in it.combinations(np.transpose(np.where(ref_img)),2):
        if (p2[0]-p1[0])**2+(p2[1]-p1[1])**2>max_sdist:
            max_sdist=(p2[0]-p1[0])**2+(p2[1]-p1[1])**2
            best_points=(p1,p2)
    deg=math.degrees(math.atan2(best_points[1][0]-best_points[0][0], best_points[1][1]-best_points[0][1]))
    deg=90-deg
    return deg%360, (deg+180)%360

def main(args):

    cgs=[]
    projs=[]
    for file_ in args.files:
        cgs.append(ftmc.CoarseGrainRNA(file_))
        try:
            projs.append(ftmp.Projection2D(cgs[-1]))
        except ValueError:
            projs.append(None)

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
        pcount=0
        for i in range(0,len(cgs)-diff):
            count+=1
            try:
                vrs1 = np.array([x for p in sorted(projs[i]._coords.keys()) for x in projs[i]._coords[p]])
                vrs2 = np.array([x for p in sorted(projs[i+diff]._coords.keys()) for x in projs[i+diff]._coords[p]])
                prmsd+=ftur.centered_rmsd(vrs1, vrs2)
                pcount+=1
            except:
                pass
            cgrmsd+=ftuc.cg_rmsd(cgs[i],cgs[i+diff])
        if count:
            cg_rmsds[diff]=cgrmsd/count
        if pcount:
            p_rmsds[diff]=prmsd/pcount
    print "projection RMSDs:", p_rmsds
    print "3D-RMSDs:", cg_rmsds
    import matplotlib.pyplot as plt
    if args.target_structure and args.hausdorff_reference:
        fig,(ax,axT, axH)=plt.subplots(3)
    elif args.target_structure:
        fig,(ax,axT)=plt.subplots(2)
    elif args.hausdorff_reference:
        fig,(ax,axH)=plt.subplots(2)
    else: 
        fig,ax=plt.subplots()

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
                try:
                    vrs1 = np.array([x for p in sorted(proj._coords.keys()) for x in proj._coords[p]])
                    prmsd=ftur.centered_rmsd(vrs1, target_vrs)
                    target_proj_rmsds.append(prmsd)
                except:
                    target_proj_rmsds.append(float("nan"))
        target_vrs = ftug.bg_virtual_residues(target)
        for cg in cgs:
            vrs1 = ftug.bg_virtual_residues(cg)
            cgrmsd=ftur.centered_rmsd(vrs1,target_vrs)
            target_rmsds.append(cgrmsd)
        xval=np.arange(len(cgs))*args.step_size
        if target_proj_rmsds:
            axT.plot(xval,target_proj_rmsds,label="Projection", color="green")
        axT.plot(xval,target_rmsds,label="3D", color="blue")
        axT.set_title("RMSD to target structure")
        axT.set_xlabel("step")
        axT.set_ylabel("RMSD")
        leg=axT.legend(framealpha=0.8, fancybox=True)
        leg.get_frame().set_linewidth(0.0)
        plt.subplots_adjust(hspace=0.4)
    if args.hausdorff_reference:
        ref_img=scipy.ndimage.imread(args.hausdorff_reference)
        ref_quarter=(scipy.ndimage.zoom(ref_img, 0.3)>150)
        degrees=get_refimg_longest_axis(ref_img)
        global_optima=[]
        local_optima=[]
        for cg in cgs:            
            score, img, params = fph.try_parameters(ref_img, args.hausdorff_scale, cg, degrees)
            local_optima.append(score)
            s, i, params = fph.globally_minimal_distance(ref_quarter, args.hausdorff_scale, cg,  
                                                          start_points=40, 
                                                          starting_rotations=degrees,
                                                          virtual_atoms=False) 
            score, img, params = fph.locally_minimal_distance(ref_img, args.hausdorff_scale, cg, 
                                                              params[1], params[2], params[0], 
                                                              maxiter=200)
            global_optima.append(score)

        axH.plot(xval,global_optima, label="Global Minimization", color="red")        
        axH.plot(xval,local_optima, label="Local Minimization", color="lavender")
        axH.set_title("Hausdorff Distance to reference image")
        axH.set_xlabel("step")
        axH.set_ylabel("Hausdorff Distance")
        leg=axH.legend(framealpha=0.8, fancybox=True)
        leg.get_frame().set_linewidth(0.0)
        plt.subplots_adjust(hspace=0.4)

    p_rmsds=np.array(p_rmsds.items())
    cg_rmsds=np.array(cg_rmsds.items())
    ax.plot(p_rmsds[:,0]*args.step_size, p_rmsds[:,1], label="Projection",color="green", marker="o")
    ax.plot(cg_rmsds[:,0]*args.step_size, cg_rmsds[:,1], label="3D", color="blue", marker="o")
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
