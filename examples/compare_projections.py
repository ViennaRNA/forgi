#!/usr/bin/python

from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
__metaclass__ = type #New style classes in Python 2.x

import argparse, random, math, sys
import itertools as it
import numpy as np
import os.path

import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.model.coarse_grain as ftmc
import forgi.projection.projection2d as ftmp
import forgi.projection.hausdorff as fph

import matplotlib.pyplot as plt
import scipy.signal
import scipy.ndimage
import scipy.misc

import multiprocessing
from functools import partial


def async_calculation((ref_img, ref_box), numFiles, filenames):
    try:
        distances=[]
        for j in range(numFiles):
            cg=ftmc.CoarseGrainRNA(filenames[j])
            distance, img, _ = fph.globally_minimal_distance(ref_img, ref_box[1]-ref_box[0], cg)
            distances.append(distance)
        return distances
    except KeyboardInterrupt:
        print("In worker {}: Keyboard Interrupt".format(id(ref_img)))
        return 

def parallel_localOpt((ref_img, ref_box), numFiles, filenames):
    try:
        distances=[]
        for j in range(numFiles):
            cg=ftmc.CoarseGrainRNA(filenames[j])
            distance, img, _ = fph.locally_minimal_distance(ref_img, ref_box[1]-ref_box[0], cg)
            distances.append(distance)
        return distances
    except KeyboardInterrupt:
        print("In worker {}: Keyboard Interrupt".format(id(ref_img)))
        return

def get_parser():
    """
    Here all commandline & help-messages arguments are defined.

    :returns: an instance of argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    #Argument(s)
    parser.add_argument('files', nargs="+", help='One or more "*.cg"/"*.coord" files.')
    #Options
    parser.add_argument('--reference', type=str, help='A *.cg file containing a "project" line or a png. Used as a reference against which all files are compared.'
                                                      'If not present, compare everything against everything. Note that only square png images are supported.')    
    parser.add_argument('--global-search', action='store_true', help='Perform global optimization over all possible projection directions. '
                                                              'Note: This is default for every cg file that does not contain a "project" entry.')
    parser.add_argument('--dpi', type=int, help='Dots per Image cloumn. Number of cells for rasterization per column and row. If reference is a png image, this is not required and would trigger rescaling of the image.')
    parser.add_argument('--scale', type=int, help='The width of the whole image in Angstrom. This is required if reference is a png image.')
    parser.add_argument('--outfile', type=str, help='Currently only used if no reference is present. Write the resulting distance matrix to this file.')
    parser.add_argument('--show', action='store_true', help='Show a figure. The figure does not show all calculated distances.')
    parser.add_argument('--cores', type=int, default=1, help='Only used for many to many comparison. Number of cores for paralellization.')
    return parser

# Parser is made available if this file is loaded as a module, to allow
# the sphinxcontrib.autoprogram extension to work properly
parser = get_parser()
if __name__=="__main__":
    args = parser.parse_args()
    if args.reference:
        if args.reference[-3:]=="png":
            ref_img=scipy.misc.imread(args.reference, flatten=True)
            if args.dpi:
                ref_img=scipy.misc.imresize(ref_img,(args.dpi+1,args.dpi+1), "nearest")
            if args.scale:
                s=args.scale/2
                ref_box=-s/2, s/2, -s/2, s/2
            else:
                parser.error("--scale is required if the reference is a png image")
            if args.show:
                fig, ax=plt.subplots()
                ax.imshow(ref_img, interpolation="none", cmap='gray')
                plt.show()
        else:
            ref_cg=ftmc.CoarseGrainRNA(args.reference)    
            try:
                ref_proj=ftmp.Projection2D(ref_cg, project_virtual_atoms=True)
            except ValueError:
                  parser.error('The reference *.cg file needs a "project" line.')
            if args.scale:            
                ref_box=fph.get_box(ref_proj, args.scale)
            else:
                ref_box=ref_proj.get_bounding_square(margin=50)
                args.scale=ref_box[1]-ref_box[0]
            if args.dpi:
                ref_img, _=ref_proj.rasterize(args.dpi, bounding_square=ref_box)    
            else:
                parser.error("If the reference is not a png image, --dpi is required")
        ref_img=(ref_img>np.zeros_like(ref_img)) #Make it a boolean array
        for f in args.files:
            cg=ftmc.CoarseGrainRNA(f)
            if args.global_search or cg.project_from is None:
                distance, img, params = fph.globally_minimal_distance(ref_img, args.scale, cg)
                fname = os.path.basename(f)
                print ("{}:\t{} distance (projected from {}, rotated by {}, offset {}. Globally optimized)".format(fname, distance, params[0], params[1], params[2] ))
            else:
                distance, img, params = fph.try_parameters(ref_img, args.scale, cg)
                fname = os.path.basename(f)
                print ("{}:\t{} distance (projected from {}, rotated by {}, offset {}. Locally optimized)".format(fname, distance, params[0], params[1], params[2]))
        if args.show:
            fig, ax=plt.subplots(2)
            ax[0].imshow(ref_img, interpolation="none", cmap='gray')
            try:
                ax[1].imshow(img, interpolation="none", cmap='gray')
            except TypeError: pass # img is probably None
            ax[0].set_title("Reference")
            ax[1].set_title("{}: {} distance".format(fname, distance))
            plt.show()
    else: #No reference given. Compare all pairs.
        if not args.dpi:
            parser.error("--dpi is required if no reference is given!")
        numFiles=len(args.files)
        distances=np.full((numFiles, numFiles), np.inf)
        combinations=numFiles*numFiles
        toProcess=[]
        for i in range(numFiles):
            cg1=ftmc.CoarseGrainRNA(args.files[i])
            try:
                ref_proj=ftmp.Projection2D(cg1, project_virtual_atoms=True)
            except ValueError:
                continue
            if args.scale:            
                ref_box=fph.get_box(ref_proj, args.scale)
            else:
                ref_box=ref_proj.get_bounding_square(margin=50)
            ref_img, _=ref_proj.rasterize(args.dpi, bounding_square=ref_box)
            ref_img=(ref_img>np.zeros_like(ref_img)) #Make it a boolean array
            #if args.show:
            #    fig, ax=plt.subplots()
            #    ax.imshow(ref_img, interpolation="none", cmap='gray')
            #    ax.set_title(args.files[i])
            #    plt.show()
            toProcess.append((ref_img, ref_box))
        if args.global_search:
            partial_calculation=partial(async_calculation, numFiles=numFiles, filenames=args.files)
        else:
            partial_calculation=partial(parallel_localOpt, numFiles=numFiles, filenames=args.files)
        pool = multiprocessing.Pool(args.cores)              
        try:
            for i, line in enumerate(pool.imap(partial_calculation, toProcess)):
                sys.stderr.write('\rCalculation {}/{} done'.format(i+1,numFiles))
                sys.stderr.flush()
                distances[i]=line

            """for j in range(numFiles):
                    a=i*numFiles+j+1
                    sys.stdout.flush()
                    sys.stdout.write("\rPerforming comparison {} of {}".format(a,combinations))
                    cg=ftmc.CoarseGrainRNA(args.files[j])
                    distance, img, _ = fph.globally_minimal_distance(ref_img, ref_box[1]-ref_box[0], cg)
                    distances[i,j]=distance"""
        except BaseException as e:
            print("Programm crashing because of a {}".format(type(e)))
            print("Distances calculated so far:")
            np.set_printoptions(threshold='nan')
            print(", ".join(os.path.basename(x) for x in args.files))
            print(distances)
            pool.terminate()
            print("terminated pool")
            raise
        finally:
            pool.close()
            pool.join()
        if args.show:
            fig, ax=plt.subplots(2)
            ax[0].imshow(ref_img, interpolation="none", cmap='gray')
            ax[1].imshow(img, interpolation="none", cmap='gray')
            ax[1].set_title("{} distance".format(distance))
            plt.show()
        print("=========================")
        np.set_printoptions(threshold='nan')
        print(distances)        
        if args.outfile:
            np.savetxt(args.outfile, distances , delimiter=",", header=", ".join(os.path.basename(x) for x in args.files))



