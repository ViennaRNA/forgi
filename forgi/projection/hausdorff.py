from __future__ import print_function
import sys, math
import numpy as np

__all__=["offsets", "modified_hausdorff_distance", "hausdorff_distance"]
__author__ = "Bernhard Thiel"
__copyright__ = "Copyright 2016"
__maintainer__ = "Bernhard Thiel"
__email__ = "thiel@tbi.univie.ac.at"


"""
A module for grid-based hausdorff distances.
This module calculates the distance between two boolean matrices.
These distances compare the location of True-values in the two matrices.
"""

def norm(offsets):
    """
    Euclidean norm.
    The following properties are used: norm(a,x)>=norm(a,0) and norm((a,0))==a
    :param offset: A tuple of length 2.
    """
    # @profile: 17% of the time spent in this function are due to the sqrt.
    #           Seems ok, because this function uses only little total runtime.
    return math.sqrt(offsets[0]**2+offsets[1]**2)


def offsets(to_iterate=[]):
    """
    An iterator over offsets (dx,dy) in the order of increasing norm((dx,dy))
    dx and dy are integers, starting at (0,0).
    """
    # We take advatage of mutable default parameter: to_iterate=[] in the function declaration.
    # If we iterate several times over offsets(), we only build to_iterate once!
    i=0
    while True:
        # Iterate over to_iterate.
        # If we reached the end, increase to_iterate and 
        # continue iteration with first added element.
        while i < len(to_iterate):
            yield to_iterate[i]            
            i+=1
        increase_range(to_iterate)

def increase_range(to_iterate):
    """
    Increase the range of offsets, which are ordered by the given norm.

    :param to_iterate: The list of offsets to be increased. 
                       A list of coordinates (dx,dy), ordered by increasing norm((dx,dy))
    """
    dists={}
    if to_iterate:
        length=norm(to_iterate[-1])
        newlength=int(length+25) #25 could be anything.
        for dx in range(-newlength, newlength+1):
            for dy in range(-newlength, newlength+1):
                dd=(dx,dy)
                if dd not in dists:
                    dists[dd]=norm(dd)
        to_iterate+=[ x for x 
                      in sorted(dists.keys(), key=lambda x: dists[x]) 
                      if length<dists[x]<newlength
                    ]
    else:
        length=25 #25 could be anything.
        for dx in range(-length, length+1):
            for dy in range(-length, length+1):
                dd=(dx,dy)
                dists[dd]=norm(dd)
        to_iterate+=[ x for x in sorted(dists.keys(), key=lambda x: dists[x]) if dists[x]<length]


def hausdorff_helperdist(p, img):
    """
    Returns the shorthest distance from a given point p to any non-zero cell in img

    :param p: a point in matrix coordinates
    :param img: A binary matrix
    """
    for (dx,dy) in offsets():
        if p[0]+dx in range(len(img)) and p[1]+dy in range(len(img[0])):
            if img[p[0]+dx, p[1]+dy]:
                return math.sqrt(dx**2+dy**2)
    else:
        return float('inf')

def hausdorff_distance(img, ref_img):
    """
    Return the grid-based Hausdorff distance between two aligned boolean matrices.
    """
    #Source: https://de.wikipedia.org/wiki/Hausdorff-Metrik
    h1=max(hausdorff_helperdist([x,y], img) for x,y in np.transpose(np.where(ref_img) ))
    h2=max(hausdorff_helperdist([x,y], ref_img) for x,y in np.transpose(np.where(img) ))
    return max( h1, h2)

def modified_hausdorff_distance( img, ref_img):
    """
    Return the grid-based Modified Hausdorff distance between two aligned boolean matrices.
    This distance was proposed in the following paper.
    It uses the mean of all distances between a po
    """
    #Source: 
    h1=sum(hausdorff_helperdist([x,y], img) 
          for x,y in np.transpose(np.where(ref_img)))/np.sum(ref_img)
    h2=sum(hausdorff_helperdist([x,y], ref_img) 
           for x,y in np.transpose(np.where(img)))/np.sum(img)
    return max( h1, h2)
