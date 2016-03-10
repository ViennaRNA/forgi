from __future__ import print_function, absolute_import
import sys, math
import numpy as np
from ..threedee.utilities import vector as ftuv
from . import projection2d as fhp
import random
import itertools as it

__all__=["offsets", "modified_hausdorff_distance", "hausdorff_distance", 
         "locally_minimal_distance", "globally_minimal_distance"]
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
    An iterator over offsets and their length ((dx,dy), norm((dx,dy))) in the order 
    of increasing norm((dx,dy))
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
                       A list of coordinates and lengths ((dx,dy), norm((dx,dy))),
                       ordered by increasing norm((dx,dy))
    """
    dists={}
    if to_iterate:
        length=norm(to_iterate[-1][0])
        newlength=int(length+25) #25 could be anything.
        for dx in range(-newlength, newlength+1):
            for dy in range(-newlength, newlength+1):
                dd=(dx,dy)
                if dd not in dists:
                    dists[dd]=norm(dd)
        to_iterate+=[ (x, dists[x]) for x 
                      in sorted(dists.keys(), key=lambda x: dists[x]) 
                      if length<dists[x]<newlength
                    ]
    else:
        length=25 #25 could be anything.
        for dx in range(-length, length+1):
            for dy in range(-length, length+1):
                dd=(dx,dy)
                dists[dd]=norm(dd)
        to_iterate+=[ (x, dists[x]) for x in sorted(dists.keys(), key=lambda x: dists[x]) if dists[x]<length]


def hausdorff_helperdist(p, img, cutoff=float("inf")):
    """
    Returns the shorthest distance from a given point p to any non-zero cell in img

    :param p: a point in matrix coordinates
    :param img: A binary matrix
    :param cutoff: If the distance is larger cutoff, return float("inf"). 
                   Increases execution speed
    """
    for dd, n in offsets():
        if n>cutoff:
            return float("inf")
        if p[0]+dd[0] in range(len(img)) and p[1]+dd[1] in range(len(img[0])):
            if img[p[0]+dd[0], p[1]+dd[1]]:
                return n
    else:
        return float('inf')


def hausdorff_distance(img, ref_img, cutoff=float("inf")):
    """
    Return the grid-based Hausdorff distance between two aligned boolean matrices.
    """
    #Source: https://de.wikipedia.org/wiki/Hausdorff-Metrik
    h1=max(hausdorff_helperdist([x,y], img, cutoff) for x,y in np.transpose(np.where(ref_img) ))
    if h1==float("inf"):
        return float("inf")
    h2=max(hausdorff_helperdist([x,y], ref_img, cutoff) for x,y in np.transpose(np.where(img) ))
    return max( h1, h2)


def hausdorff_distance_old(img, ref_img, cutoff=float("inf")):
    """
    Return the grid-based Hausdorff distance between two aligned boolean matrices.
    """
    #Source: https://de.wikipedia.org/wiki/Hausdorff-Metrik
    h1=max(hausdorff_helperdist([x,y], img, cutoff) for x,y in np.transpose(np.where(ref_img) ))
    h2=max(hausdorff_helperdist([x,y], ref_img, cutoff) for x,y in np.transpose(np.where(img) ))
    return max( h1, h2)

def modified_hausdorff_distance( img, ref_img, _=None):
    """
    Return the grid-based Modified Hausdorff distance between two aligned boolean matrices.
    This distance was proposed in the following paper.
    It uses the mean of all distances instead of the max.
    """
    #Source: 
    h1=sum(hausdorff_helperdist([x,y], img) 
          for x,y in np.transpose(np.where(ref_img)))/np.sum(ref_img)
    h2=sum(hausdorff_helperdist([x,y], ref_img) 
           for x,y in np.transpose(np.where(img)))/np.sum(img)
    return max( h1, h2)

def to_polar(x):
  return np.array(ftuv.spherical_cartesian_to_polar(x))

def from_polar(x):
  return ftuv.spherical_polar_to_cartesian(x)

def get_box(projection, width, offset=np.array((0,0))):
    left, right, down, up=projection.get_bounding_square()
    center_hor=left+(right-left)/2
    center_ver=down+(up-down)/2
    box=(center_hor-width/2+offset[0], center_hor+width/2+offset[0], 
        center_ver-width/2+offset[1], center_ver+width/2+offset[1])
    return box

def compare(ref_img, cg, scale, projection_direction, distance, 
            rotation, offset, cutoff=float("inf"), nfev=[0], write_fev=False):
    """
    A helper function for locally_minimal_distance_one.
    For parameter: See documentation of locally_minimal_distance
    :param projection_direction: in polar coordinates.
    :param rotation: The in plane rotation of the projection in degrees
    :param offset: The translation of the proj. image to the reference image in Angstrom (dx,dy)
    :param nfev: A list. The last list element is incremented every time this function 
                 is called to keep track of the number of function evaluations.
    :param write_fev: If true, also print the number of the last element of nfev.
    """
    nfev[-1]+=1
    if write_fev:
        print("Number of function evaluations: ", nfev[-1])
        nfev[-1]=0
    proj=fhp.Projection2D(cg, from_polar([1]+list(projection_direction)))
    proj.rotate(rotation)
    box=get_box(proj, scale, offset)
    img,_=proj.rasterize(len(ref_img)-1, bounding_square=box, warn=False)
    # Convert img to a boolean image
    img=(img>np.zeros_like(img))
    return distance(ref_img, img, cutoff), img
    
def locally_minimal_distance_one(ref_img, scale, cg, start_rot, proj_dir, maxiter, distance):
    """
    A helper function for locally_minimal_distance.
    For parameters see documentation of locally_minimal_distance
    :startrot: The starting rotation in the projection plane (in degrees).
    """
    curr_best_rotation=start_rot
    curr_best_offs=np.array((0,0))
    if proj_dir is not None:
        curr_best_pro=proj_dir
    else:
        curr_best_pro=to_polar(cg.project_from)[1:]        
    curr_best_score, _ = compare(ref_img, cg, scale, curr_best_pro, distance=distance,
                                 rotation=curr_best_rotation, offset=curr_best_offs) 
    #Stepwise improve score
    for i in range(maxiter):
        found=False
        # ---------------
        # Optimize offset      
        # --------------- 
        change_offs=np.array((0,0)) #Angstrom dx,dy
        for co in [np.array((1,0)), np.array((0,1)), np.array((-1,0)), np.array((0,-1))]:
            for i in range(10):
                tmp_score, _ = compare(ref_img, cg, scale, curr_best_pro, distance=distance,
                                       rotation=curr_best_rotation, offset=curr_best_offs+co,
                                       cutoff=curr_best_score)
                if tmp_score>curr_best_score:
                    break
                elif tmp_score<curr_best_score:
                    curr_best_score=tmp_score
                    change_offs=co
                    break
                else:
                   #Function value did not change. Our stepwidth was to small.
                   co=co*2
            else:
                pass #print ("Offset: No change in direction ", co)
        if np.all(change_offs==np.array((0,0))):
            found=True
        #Save the new currently best offset
        curr_best_offs=curr_best_offs+change_offs

        # -----------------
        # Optimize rotation      
        # -----------------
        change_rot=0 #degrees
        for cr in (-0.5, 0.5):
            for i in range(10):
                tmp_score, _ = compare(ref_img, cg, scale, curr_best_pro, distance=distance,
                                       rotation=curr_best_rotation+cr, offset=curr_best_offs,
                                       cutoff=curr_best_score)
                if tmp_score>curr_best_score:
                    break
                elif tmp_score<curr_best_score:
                    curr_best_score=tmp_score
                    change_rot=cr
                    break
                else:
                   #Function value did not change. Our stepwidth was to small.
                   cr=cr*2
            else:
                pass #print ("Rotation: No change in direction ", cr)
        if change_rot!=0:
            # If we didn't change the offset and didn't change the rotation, found stays True
            found=False
        curr_best_rotation=curr_best_rotation+change_rot

        # -----------------------------
        # Optimize projection direction      
        # -----------------------------
        change_pro=np.array((0.,0.)) #polar coordinates
        for cp in [(0,0.002), (0.002,0), (0,-0.002), (-0.002,0),(0.002,0.002),
                                  (-0.002,-0.002), (0.002,-0.002),(-0.002,0.002)]:
            for i in range(10):
                tmp_score, _ = compare(ref_img, cg, scale, curr_best_pro+cp, distance=distance, 
                                       rotation=curr_best_rotation, offset=curr_best_offs,
                                       cutoff=curr_best_score)
                if tmp_score>curr_best_score:
                    break
                elif tmp_score<curr_best_score:
                    curr_best_score=tmp_score
                    change_pro=cp
                    break
                else:
                   #Function value did not change. Our stepwidth was to small.
                   cp=np.array(cp)*2
            else:
                pass #print ("Projection: No change in direction ", cp)
        if not np.all(change_pro==np.array((0.,0.))):
            # If we didn't change offset, projection or rotation, found still stays True
            found=False
        curr_best_pro=curr_best_pro+change_pro
        if found:
            break
    else:
        pass #print("Max-iter reached", maxiter)
    best_score, img = compare(ref_img, cg, scale, curr_best_pro, distance=distance,
                              rotation=curr_best_rotation, offset=curr_best_offs,
                              cutoff=curr_best_score)#, write_fev=True)
    assert best_score==curr_best_score
    return curr_best_score, img, curr_best_pro

def locally_minimal_distance(ref_img, scale, cg, proj_dir=None, maxiter=50, 
                             distance=hausdorff_distance):
    """
    :param ref_img: The reference image. A boolean square matrix.
    :param scale: The edge length in Angstrom of the reference image.
    :param cg: The coarse grain RNA to match to the projection.
    :param proj_dir: The starting projection direction in spherical polar coordinates.
                     If this is not given: uses cg.project_from.
    :param distance: a function with signature like hausdorff_distance 
                     or modified_hausdorff_distance.
    :returns: A triple: score, optimal_img, projection_direction
    """
    best_score=float("inf")
    #Estimate the best in-plane rotation (without local optimization)        
    if proj_dir is  None:
        proj_dir=to_polar(cg.project_from)[1:]        
    for start_rot in [0, 180]:

        initial_score, _ = compare(ref_img, cg, scale, proj_dir, distance=distance,
                                   rotation=start_rot, offset=np.array((0,0)))       
        if initial_score<best_score:
            best_score=initial_score
            best_rot=start_rot
    #Now minimize locally
    best_score, best_img, best_pro=locally_minimal_distance_one(ref_img, scale, cg, best_rot, 
                                                     proj_dir, int(maxiter/2+1), distance)   
    return best_score, best_img, best_pro


def get_start_points(numPoints):
    """
    Return numPoints equally-distributed points on the unit sphere in polar coordinates.

    Implements the 2nd algorithm from https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
    """
    a=4*math.pi/numPoints
    d=math.sqrt(a)
    Mt=int(math.pi/d)
    dt=math.pi/Mt
    df=a/dt
    points=[]
    for m in range(Mt):
        theta=math.pi*(m+0.5)/Mt
        Mf=int(2*math.pi*math.sin(theta)/df)
        for n in range(Mf):
            phi=2*math.pi*n/Mf
            points.append(np.array((theta, phi)))
    return points

def get_longest_img_diameter(img, scale):
    maxl=0
    for x1,y1 in np.transpose(np.where(img)):
        for x2,y2 in np.transpose(np.where(img)):
            l=(x1-x2+1)**2+(y1-y2+1)**2
            if l>maxl:
                maxl=l
    return math.sqrt(maxl)*scale/len(img)

def globally_minimal_distance(ref_img, scale, cg, start_points=20, maxiter=5, 
                             distance=hausdorff_distance):
    """
    See parameters of locally_minimal_distance.
    Calls locally_minimal_distance for several starting projection directions.
    Uses several Heuristics to speed up the process.
    :param maxiter: Maximal iteration for each local optimization
    :param start_points: Number of starting projection directions. For each, local
                         optimization with maxiter steps is performed.
    """
    best_score=float("inf")
    decrease=float("inf")
    sp=get_start_points(start_points)

    #We want to cover many different directions with the frist few iterations.
    random.seed(1)
    random.shuffle(sp)
    #Longest extention in the image
    longest_distance_image=get_longest_img_diameter(ref_img, scale)
    if longest_distance_image==0:
        raise ValueError("Reference image not working. Probably empty")
    la_heur=0
    no_heur=0
    score_heur=0
    for project_dir in sp:
        proj=fhp.Projection2D(cg, from_polar([1]+list(project_dir)))
        if abs(proj.get_longest_leaf_leaf_distance()-longest_distance_image)>4*scale/len(ref_img): #4 pixels is arbitrary heuristic
            la_heur+=1
            continue
        initial_score, _ = compare(ref_img, cg, scale, project_dir, distance=distance,
                                   rotation=0, offset=np.array((0,0)))        
        #print("Score ", initial_score)
        if initial_score>best_score+decrease*1.5:
            #print("cnt (cutoff = ",best_score+decrease*1.5 ," )")
            score_heur+=1
            continue
        no_heur+=1
        score, img, pro=locally_minimal_distance( ref_img, scale, cg, project_dir, 
                                                  maxiter, distance )
        #print ("optimized to ", score, " (curr_best is {})".format(best_score))
        if score<best_score:
            best_score=score
            best_img=img
            best_pro=pro
            decrease=initial_score-score
    if no_heur==0:
        return float("inf"), None, None
    score, img, pro=locally_minimal_distance( ref_img, scale, cg, best_pro, 
                                                  20*maxiter, distance )
    #print("Calculations performed: {}, skipped (longest distance): {}, skipped (score): {}".format(no_heur, la_heur, score_heur))
    return score, img, pro



