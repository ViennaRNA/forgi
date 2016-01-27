from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

import sys, random, math
import itertools as it
import collections as col

import numpy as np
import matplotlib.pyplot as plt

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.projection2d as ftmp
import forgi.threedee.utilities.vector as ftuv

import argparse

def get_parser():
    """
    Here all commandline & help-messages arguments are defined.

    :returns: an instance of argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    #Argument(s)
    parser.add_argument('cgfiles', nargs='+', help='One or more *.cg/*.coord files holding the RNA to plot.')
    parser.add_argument('--show-direction', type=bool, default=False, help='Print the projection direction in the plot')
    parser.add_argument('--out', '-o', type=str, nargs='+', help='One or more outfiles to save the resulting figure to. '
                        'The file format will be determined by the file ending. Supported formats are "svg", "png" and "pgf"')
    parser.add_argument('--out-path', type=str, nargs='?', help='Optional path, used for all files given with the "--out" option')
    parser.add_argument('--show-filename', type=bool, default=False, help='Print the filename of the input file in the figure')




#The CONDENSE constant is used to simulate a reduced resolution. 
#The algorithm for doing so is work in progress.
#For now, a value less than 20 is reasonable.
CONDENSE=0

if __name__=="__main__":
    #This script expects a list of *.cg/ *.coord files as command line arguments (and nothing else)
    files=sys.argv[1:]

    #Uncomment the following line to display the files in a random order.
    #random.shuffle(files)

    #Prepare the pyplot figure
    totalFigures=len(files)
    figuresPerLine=int(math.ceil(math.sqrt(totalFigures)))
    fig, ax=plt.subplots(int(math.ceil(totalFigures/figuresPerLine)),figuresPerLine, squeeze=False)
    
    #Uncomment the following line to change the background color of the figure (not the plot).
    #fig.patch.set_facecolor('black')

    #Plot one projection per file.
    for i, file_ in enumerate(files):
        #get the subplot axes (Note: axes != axis in matplotlib)
        current_axes=ax[i//figuresPerLine, i%figuresPerLine]

        #Parse the file
        cg=ftmc.CoarseGrainRNA(file_)

        # Random projection direction. Change to direction=[1.,1.,0.] to set a specific direction
        direction= [-0.69562761, -0.19249465, -0.69213296]#ftuv.get_random_vector()

        #Generate the projection object
        proj=ftmp.Projection2D(cg, direction, rotation=180)   

        #Simulate a reduced resolution of the image.     
        proj.condense(CONDENSE)
        
        elems=list(proj._coords.keys())
        random.shuffle(elems)
        hairpins=[ x for x in elems if x[0]=="h" ]
        multiloops=[ x for x in elems if x[0]=="m"]
        elems=hairpins[:2]+[multiloops[0]]+[elems[0]]
        elems=["h1", "h8", "m31", "h6", "h5", "s44", "m12", "m15", "h16"]
        comb=list(it.combinations(elems, 2))
        #Plot the projection #
        proj.plot(ax[i//figuresPerLine, i%figuresPerLine], margin=15, linewidth=5, add_labels=set(elems), line2dproperties={"color":"darkgray", "linestyle":"-"},
                  show_distances=comb, print_distances=True)

        #Uncomment to set a substring of the filename as a title
        #current_axes.set_title(file[-15:])

        #Hide the x- and y axis.
        current_axes.get_xaxis().set_visible(False)
        current_axes.get_yaxis().set_visible(False)

        #Uncomment the following lines to print the projection direction and the filename in the plot.
        current_axes.text(0.01,0.01,"Projection direction: ({},{},{})".format(round(direction[0],3), round(direction[1],3), round(direction[2],3)), transform=current_axes.transAxes)
        current_axes.text(0.01,0.97,"File: {}".format(file_), transform=current_axes.transAxes)
        #Uncomment the following line to change the backgroundcolor of the plot area.
        #current_axes.set_axis_bgcolor('black')

    
    #Hide additional subplots with no projection on them.
    for i in range(len(files),int(math.ceil(totalFigures/figuresPerLine))*figuresPerLine):
        ax[i//figuresPerLine, i%figuresPerLine].axis('off')

    # Reduce the space outside of the plots and between the subplots.
    plt.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.975, wspace=0.05, hspace=0.05)
    
    #Show the plot and clear it from the internal memory of matplotlib.
    plt.show()
