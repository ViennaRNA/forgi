from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

import sys, random, math, os
import itertools as it
import collections as col

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import forgi.threedee.model.coarse_grain as ftmc
import forgi.projection.projection2d as ftmp
import forgi.threedee.utilities.vector as ftuv

import argparse

#################### MAIN FUNCTION #####################################

def main(args):    
    files=args.cgfiles

    #Uncomment the following line to display the files in a random order.
    #random.shuffle(files)

    #Prepare the pyplot figure
    totalFigures=len(files)
    figuresPerLine=int(math.ceil(math.sqrt(totalFigures)))
    fig, ax=plt.subplots(int(math.ceil(totalFigures/figuresPerLine)),figuresPerLine, squeeze=False, figsize=(8,8))
    
    #Background color of figure (not plot)
    if args.style=="WOB":
        fig.patch.set_facecolor('black')

    #Plot one projection per file.
    for i, file_ in enumerate(files):
        #get the subplot axes (Note: axes != axis in matplotlib)
        current_axes=ax[i//figuresPerLine, i%figuresPerLine]

        #Parse the file
        cg=ftmc.CoarseGrainRNA(file_)

        # Random projection direction, if no direction present in the file
        if args.proj_direction:
            direction=list(map(float, args.proj_direction.split(",")))
        elif cg.project_from is not None:
            direction=cg.project_from
        else:
            direction=ftuv.get_random_vector()

        #Generate the projection object
        proj=ftmp.Projection2D(cg, direction, rotation=180, project_virtual_atoms=args.virtual_atoms)   

        #Simulate a reduced resolution of the image.
        if args.condense:
            proj.condense(args.condense)

        target_elems=[]
        if args.show_distances:                

            try:
                num_elems=int(args.show_distances)
            except:
                target_elems=args.show_distances.split(",")
            else:
                if num_elems>len(proj._coords.keys()):
                    raise ValueError("--show-distances must not be greater {} for the current projection ({}:'{}')".format(len(proj._coords.keys()), i, file_))
                elems=list(proj._coords.keys())
                random.shuffle(elems)
                while len(target_elems)<num_elems:
                    r=random.random()
                    if r<0.4:
                        hairpins=[ x for x in elems if x[0]=="h" and x not in target_elems ]
                        if hairpins:
                            target_elems.append(hairpins[0])
                            continue
                    if r<0.6:                        
                        multiloops=[ x for x in elems if x[0]=="m" and x not in target_elems ]
                        if multiloops:
                            target_elems.append(multiloops[0])
                            continue
                    others=[ x for x in elems if x not in target_elems]
                    target_elems.append(others[0])
        comb=list(it.combinations(target_elems, 2))
        #print(comb, target_elems)
        if args.label_elements:
            target_elems=list(proj._coords.keys())
        line2dproperties={}
        if args.style=="BOW":
            line2dproperties["color"]="black"
        elif args.style=="WOB":
            line2dproperties["color"]="white"

        #Plot the projection #
        proj.plot(current_axes, margin=15, linewidth=3, add_labels=set(target_elems), line2dproperties=line2dproperties,
                  show_distances=comb, print_distances=args.print_distances)

        #Uncomment to set a substring of the filename as a title
        #current_axes.set_title(file[-15:])

        #Hide the x- and y axis.
        current_axes.get_xaxis().set_visible(False)
        current_axes.get_yaxis().set_visible(False)

        #Print the projection direction and the filename in the plot.
        if args.show_direction or args.p:
            current_axes.text(0.01,0.01,"Projection direction: ({},{},{})".format(round(direction[0],3), round(direction[1],3), round(direction[2],3)), transform=current_axes.transAxes)
        if args.show_filename or args.p:
            current_axes.text(0.01,0.99,"File: {}".format(file_), transform=current_axes.transAxes, verticalalignment='top',)

        #Change the backgroundcolor of the plot area.
        if args.style=="WOB":
            current_axes.set_axis_bgcolor('black')

    
    #Hide additional subplots with no projection on them.
    for i in range(len(files),int(math.ceil(totalFigures/figuresPerLine))*figuresPerLine):
        ax[i//figuresPerLine, i%figuresPerLine].axis('off')

    # Reduce the space outside of the plots and between the subplots.
    plt.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.975, wspace=0.05, hspace=0.05)
    

    if args.out:
        for ofname in args.out:
            if args.out_path:
                ofname=os.path.join(args.out_path, ofname)
            ofname=os.path.expanduser(ofname)
            plt.savefig(ofname,format=ofname[-3:])
    if not args.out or args.show:
        #Show the plot and clear it from the internal memory of matplotlib.
        plt.show()

########################### END OF MAIN FUNCTION #########################################

def get_parser():
    """
    Here all commandline & help-messages arguments are defined.

    :returns: an instance of argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    #Argument(s)
    parser.add_argument('cgfiles', nargs='+', help='One or more *.cg/*.coord files holding the RNA to plot.')
    parser.add_argument('--show-direction', action="store_true",help='Print the projection direction in the plot')
    parser.add_argument('--proj-direction', type=str, help='Use the given projection direction instead of the one from the file. A comma seperated triple of floats (with no whitespace)')      
    parser.add_argument('--show-filename', action="store_true", default=False, help='Print the filename of the input file in the figure')  
    parser.add_argument('-p', action="store_true", help='Shortcut for --show-direction and -- show_filename. Note that texts are not visible in all styles.')    
    parser.add_argument('--label-elements', default=False, action="store_true", help='Label all coarse-grained elements in the plot.')
    parser.add_argument('--print-distances', default=False, action="store_true", help='Print distances for all elements given in --show-distances at the side in the plot')
    parser.add_argument('--out', '-o', type=str, nargs='+', help='One or more outfiles to save the resulting figure to. '
                        'The file format will be determined by the file ending. Formats could be e.g. "svg", "png" or "pgf"')
    parser.add_argument('--out-path', type=str, nargs='?', help='Optional path, used for all files given with the "--out" option')

    parser.add_argument('--style', type=str, default="DEF", choices=["DEF","WOB", "BOW", "COL"], 
                        help='Plot style. "DEF" (default: color on white), "WOB" (white on black), "BOW" (Black on white), "COL" same as "DEF" ')
    parser.add_argument('--condense', type=float, help='Simulate resolution reduction. This is an experimental featyure and does not work very well.')
    parser.add_argument('--show-distances', type=str, help='Either an Integer, or a ","-separated list of coarse grained elements.')
    parser.add_argument('--virtual-atoms', action="store_true", help='Show virtual atoms (Slow).')
    parser.add_argument('--show', action="store_true", help='Show the plot. If args.out is not given, this is implicitely set to true..')
    return parser
parser=get_parser()
if __name__=="__main__":
    main(parser.parse_args())
