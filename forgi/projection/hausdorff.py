from __future__ import print_function, absolute_import, division
from builtins import range
import sys
import math
import numpy as np
from ..threedee.utilities import vector as ftuv
from . import projection2d as fhp
import random
import itertools as it

__all__ = ["offsets", "modified_hausdorff_distance", "hausdorff_distance",
           "locally_minimal_distance", "globally_minimal_distance", "get_box",
           "get_longest_img_diameter", "try_parameters"]
__author__ = "Bernhard Thiel"
__copyright__ = "Copyright 2016"
__maintainer__ = "Bernhard Thiel"
__email__ = "thiel@tbi.univie.ac.at"


"""
A module for grid-based hausdorff distances.
This module calculates the distance between two boolean matrices and
the distance between a forgi.projection.projection2d.Projection2D object 
and a boolean matrix
These distances are based on the location of True-values in the two matrices.
"""

try:
    profile  # The @profile decorator from line_profiler (kernprof)
except NameError:
    def profile(x):
        return x

##############################################################################
# Helper functions for iterating the grid in "spirals"
##############################################################################


def norm(offset):
    """
    Euclidean norm.
    :param offset: A tuple of length 2.
    :returns: A float
    """
    return math.sqrt(offset[0]**2 + offset[1]**2)


def offsets(skip=0, to_iterate=[], toskip={}):  # pylint: disable=W0102
    """
    An iterator over offsets and their length ((dx,dy), norm((dx,dy))) in the order 
    of increasing norm((dx,dy))
    dx and dy are integers, starting at (0,0).

    :param skip: A float. The iterator may skip cell up to SKIP away from (0,0) 
                 before returning the first value.
    :param to_iterate: Do not use this. The mutable default parameter is used as a hack for caching.
    :param toskip: Do not use this. The mutable default parameter is used as a hack for caching.

    :yields: A tuple ((dx,dy), n) where dx and dy are integers and n=norm((dx,dy)).
    """
    # We take advatage of mutable default parameter: to_iterate=[] in the function declaration.
    # If we iterate several times over offsets(), we only build to_iterate once!
    if skip == 0:
        i = 0
    else:
        i = toskip[skip]
    while True:
        # Iterate over to_iterate.
        # If we reached the end, increase to_iterate and
        # continue iteration with first added element.
        len_to_iterate = len(to_iterate)
        while i < len_to_iterate:
            yield to_iterate[i]
            i += 1
        increase_range(to_iterate, toskip)


def increase_range(to_iterate, to_skip):
    """
    Increase the range of offsets, which are ordered by their norm. Used by offsets().

    :param to_iterate: The list of offsets to be increased. 
                       A list of coordinates and lengths ((dx,dy), norm((dx,dy))),
                       ordered by increasing norm((dx,dy))
    """
    dists = {}
    if to_iterate:
        length = to_iterate[-1][1]
        newlength = int(length + 25)  # 25 could be anything.
        for dx in range(-newlength, newlength + 1):
            for dy in range(-newlength, newlength + 1):
                dd = (dx, dy)
                if dd not in dists:
                    dists[dd] = norm(dd)
        oldlen = to_iterate[-1][0][0]**2 + to_iterate[-1][0][1]**2
        #print ("sorting {} elements for newlength {}".format(len(dists.keys()), newlength))
        for x in sorted(dists.keys(), key=lambda x: dists[x]):
            if length < dists[x] < newlength:
                to_iterate.append((x, dists[x]))
                if x[0]**2 + x[1]**2 > oldlen:
                    to_skip[x[0]**2 + x[1]**2] = len(to_iterate) - 2
    else:
        length = 25  # 25 could be anything.
        for dx in range(-length, length + 1):
            for dy in range(-length, length + 1):
                dd = (dx, dy)
                dists[dd] = norm(dd)
        oldlen = 0
        for x in sorted(dists.keys(), key=lambda x: dists[x]):
            if not dists[x] < length:
                continue
            to_iterate.append((x, dists[x]))
            if x[0]**2 + x[1]**2 > oldlen:
                to_skip[x[0]**2 + x[1]**2] = len(to_iterate) - 2

##############################################################################
# Grid based distances
##############################################################################


@profile
def hausdorff_helperdist(p, img, cutoff=float("inf"), skip=0):
    """
    Returns the shorthest distance from a given point p to any non-zero cell in img.
    This is used by hausdorff_distance()

    :param p: a point in matrix coordinates
    :param img: A binary matrix
    :param cutoff: If the distance is larger cutoff, return float("inf"). 
                   Increases execution speed
    :param skip: The function does not have to search for a distance smaller than skip 
                 (used for speedup)
    """
    dpi = len(img)
    p0 = p[0]
    p1 = p[1]
    for dd, n in offsets(skip):
        x = p0 + dd[0]
        y = p1 + dd[1]
        if 0 <= x < dpi and 0 <= y < dpi and img[x, y]:
            return n
        if n > cutoff:
            return float("inf")
    return float('inf')


def hausdorff_distance(img, ref_img, cutoff=float("inf")):
    """
    Return the grid-based Hausdorff distance between two aligned boolean matrices.

    :param img, ref_img: Two aligned boolean 2D matrices with the same shape.
    :param cutoff: A float. If the distance is greater than cutoff, return float("inf").
                   (Used to increase execution speed in certain cases)
    """
    return hausdorff_distance_new(img, ref_img, cutoff)
    # Source: https://de.wikipedia.org/wiki/Hausdorff-Metrik
    # h1=max(hausdorff_helperdist([x,y], img, cutoff) for x,y in np.transpose(np.where(ref_img) ))
    # if h1==float("inf"):
    #     return float("inf")
    # h2=max(hausdorff_helperdist([x,y], ref_img, cutoff) for x,y in np.transpose(np.where(img) ))
    # return max( h1, h2)


@profile
def hausdorff_distance_new(img, ref_img, cutoff=float("inf")):
    """
    A faster implementation using skip, called by hausdorff_distance().
    This implementation is faster for high Hausdorff distances.
    """
    x = 0
    y = 0
    oldx = None
    oldy = None
    skip = 0
    try:
        maxh = 0
        oldh = 0
        for x, y in np.transpose(np.where(ref_img)):
            if oldh == 0:
                oldh = hausdorff_helperdist([x, y], img, cutoff)
                oldx = x
                oldy = y
            else:
                # skip always has to be strictly smaller than oldh-norm(dx,dy)!
                skip = max(0, int(oldh - norm([x - oldx, y - oldy]) - 0.1))
                oldh = hausdorff_helperdist([x, y], img, cutoff, int(skip**2))
                oldx = x
                oldy = y
            if oldh > maxh:
                maxh = oldh
            if oldh == float("inf"):
                return float("inf")
        h1 = maxh
        maxh = 0
        oldh = 0
        for x, y in np.transpose(np.where(img)):
            if oldh == 0:
                oldh = hausdorff_helperdist([x, y], ref_img, cutoff)
                oldx = x
                oldy = y
            else:
                skip = max(0, int(oldh - norm([x - oldx, y - oldy]) - 0.1))
                oldh = hausdorff_helperdist(
                    [x, y], ref_img, cutoff, int(skip**2))
                oldx = x
                oldy = y
            if oldh > maxh:
                maxh = oldh
            if oldh == float("inf"):
                return float("inf")
        h2 = maxh
        return max(h1, h2)
    except KeyboardInterrupt:
        print("x: {}, y: {}, oldx: {}, oldy: {}, oldh: {}, maxh: {}, "
              "h1: {}, skip={}".format(x, y, oldx, oldy, oldh, maxh, h1, skip))
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2)
        ax[0].imshow(ref_img, interpolation="none", cmap='gray')
        ax[1].imshow(img, interpolation="none", cmap='gray')
        ax[0].set_title("Reference")
        ax[1].set_title("Image")
        plt.show()
        raise


def modified_hausdorff_distance(img, ref_img, _=None):
    """
    Return the grid-based Modified Hausdorff distance between two aligned boolean matrices.
    This distance was proposed in the following paper: TODO
    It uses the mean of all distances instead of the max.
    The current implementation is not optimized
    """
    # Source:
    h1 = sum(hausdorff_helperdist([x, y], img)
             for x, y in np.transpose(np.where(ref_img))) / np.sum(ref_img)
    h2 = sum(hausdorff_helperdist([x, y], ref_img)
             for x, y in np.transpose(np.where(img))) / np.sum(img)
    return max(h1, h2)


def tp_fp_distance(img, ref_img, _=None):
    tp = np.sum(np.logical_and(img, ref_img))
    alle = np.sum(np.logical_or(img, ref_img))
    return alle / tp


def combined_distance(img, ref_img, _=None):
    return tp_fp_distance(img, ref_img) + hausdorff_distance(img, ref_img)
##############################################################################
# Helper functions for working with projections
##############################################################################


def to_polar(x):
    return np.array(ftuv.spherical_cartesian_to_polar(x))


def from_polar(x):
    return ftuv.spherical_polar_to_cartesian(x)


def get_box(projection, width, offset=np.array((0, 0))):
    left, right, down, up = projection.get_bounding_square()
    center_hor = left + (right - left) / 2
    center_ver = down + (up - down) / 2
    box = (center_hor - width / 2 + offset[0], center_hor + width / 2 + offset[0],
           center_ver - width / 2 + offset[1], center_ver + width / 2 + offset[1])
    return box


@profile
def locally_minimal_distance(ref_img, scale, cg,
                             start_rot=None, offset=None, proj_dir=None,
                             maxiter=50, advanced=False, virtual_atoms=True,
                             distance=hausdorff_distance):
    """
    Local optimization of the Hausdorff distance between a image and a CoarseGrainRNA

    For a given reference image and a given CoarseGrainRNA, the distance is a function of
    continuouse variables (projection direction, offset and rotation), but is a stepfunction, 
    because it is grid-based. Therefor traditional optimization approaches, 
    which rely on a gradient, fail.
    This implementation tests, whether the distance is increasing or decreasing in 
    different directions, no matter how far one has to walk in that direction before 
    the distance changes. It then follows the deepest decrease (not the steepest decrease!) 

    :param ref_img: The reference image. A boolean square matrix.
    :param scale: The edge length in Angstrom of the reference image.
    :param cg: The coarse grain RNA to match to the projection.
    :param start_rot: The in-plane rotation of the projection.
    :param offset: The offset that will be applied to the projection in Angstrom. A np.array([x,y]).
    :param proj_dir: The starting projection direction in spherical polar coordinates.
                     If this is not given: uses cg.project_from.
    :param maxiter: Maximal number of iterations. An int.
    :param distance: a function with signature like hausdorff_distance 
    :param advanced: Try steps in more directions (takes longer)
    :param virtual_atoms: Whether or not to project the virtual atoms. A boolean.

    :returns: A triple: (distance, image, parameters)
              Where distance is a float,image a matrix and params is a triple:
              np.array([theta, phi]), degrees, np.array([x_offset, y_offset])
    """
    dpi = len(ref_img)
    cell_length = scale / dpi
    # Heuristical parameters
    #: The minimal stepwidth in Angstrom for offsets is determined by the width a single pixel.
    offset_stepwidth = max(1, int(cell_length / 2.5))
    #: The change in projection direction also depends on the resolution, but the correlation is less straight forward.
    scale_projectionsteps = max(1, int(cell_length / 5) * 4)
    #print("Projection Step is {}, scale {}".format(0.002*scale_projectionsteps, scale_projectionsteps))
    if start_rot:
        curr_best_rotation = start_rot
    else:
        curr_best_rotation = 0
    if offset is None:
        curr_best_offs = np.array([0, 0])
    else:
        curr_best_offs = offset
    if proj_dir is not None:
        curr_best_pro = np.array(proj_dir)
    else:
        curr_best_pro = to_polar(cg.project_from)[1:]

    # Initial projection object and initial score
    projection = fhp.Projection2D(cg, from_polar([1] + list(curr_best_pro)),
                                  project_virtual_atoms=virtual_atoms)
    box = get_box(projection, scale, curr_best_offs)
    img, _ = projection.rasterize(dpi, bounding_square=box, warn=False,
                                  rotate=curr_best_rotation)
    curr_best_score = distance(ref_img, img)

    # Stepwise improve score
    for i in range(maxiter):
        found = False
        # ---------------
        # Optimize offset
        # ---------------
        best_change_offs = np.array((0, 0))  # Angstrom dx,dy
        directions = [np.array((1, 0)), np.array(
            (0, 1)), np.array((-1, 0)), np.array((0, -1))]
        if advanced:
            directions += [np.array((1, 1)), np.array((1, -1)),
                           np.array((-1, 1)), np.array((-1, -1))]
        for co in directions:
            co = co * offset_stepwidth
            for i in range(10):
                box = get_box(projection, scale, curr_best_offs + co)
                img, _ = projection.rasterize(dpi, bounding_square=box, warn=False,
                                              rotate=curr_best_rotation)
                tmp_score = distance(ref_img, img, curr_best_score)
                if tmp_score > curr_best_score:
                    # print("co was {}".format(co)) #This was used to find the optimal stepwidth
                    break
                elif tmp_score < curr_best_score:
                    curr_best_score = tmp_score
                    best_change_offs = co
                    #print("co was {}".format(co))
                    break
                else:
                    # Function value did not change. Our stepwidth was to small.
                    co = co * 2
            else:
                pass  # print ("Offset: No change in direction ", co)
        if np.all(best_change_offs == np.array((0, 0))):
            found = True
        # Save the new currently best offset
        curr_best_offs = curr_best_offs + best_change_offs

        # The Box to be used in the following steps
        box = get_box(projection, scale, curr_best_offs)

        if not curr_best_score:
            break  # We reached zero
        # ----------------
        # Optimize rotation
        # -----------------
        best_change_rot = 0  # degrees
        for cr in (-0.5, 0.5):
            for i in range(10):
                img, _ = projection.rasterize(dpi, bounding_square=box, warn=False,
                                              rotate=curr_best_rotation + cr)
                tmp_score = distance(ref_img, img, curr_best_score)
                if tmp_score > curr_best_score:
                    #print("cr was {}".format(cr))
                    break
                elif tmp_score < curr_best_score:
                    curr_best_score = tmp_score
                    best_change_rot = cr
                    #print("cr was {}".format(cr))
                    break
                else:
                    # Function value did not change. Our stepwidth was to small.
                    cr = cr * 2
            else:
                pass  # print ("Rotation: No change in direction ", cr)
        if best_change_rot != 0:
            found = False
        else:
            pass  # If we didn't change the offset and didn't change the rotation, found stays True
        curr_best_rotation = curr_best_rotation + best_change_rot
        if not curr_best_score:
            break  # We reached zero
        # -----------------------------
        # Optimize projection direction
        # -----------------------------
        change_pro = np.array((0., 0.))  # polar coordinates
        directions = [(0, 0.002), (0.002, 0), (0, -0.002), (-0.002, 0), (0.002, 0.002),
                      (-0.002, -0.002), (0.002, -0.002), (-0.002, 0.002)]
        if advanced:
            directions += [(0.001, 0.002), (-0.001, 0.002), (0.002, 0.001), (0.002, -0.001),
                           (0.001, -0.002), (-0.001, -0.002), (-0.002, 0.001), (-0.002, -0.001)]
        for cp in directions:
            cp = (cp[0] * scale_projectionsteps, cp[1] * scale_projectionsteps)
            # print("=======================================")
            for i in range(10):
                tmp_projection = fhp.Projection2D(cg, from_polar([1] + list(curr_best_pro + cp)),
                                                  project_virtual_atoms=virtual_atoms)
                box = get_box(tmp_projection, scale, curr_best_offs)
                img, _ = tmp_projection.rasterize(dpi, bounding_square=box, warn=False,
                                                  rotate=curr_best_rotation)
                tmp_score = distance(ref_img, img, curr_best_score)
                if tmp_score > curr_best_score:
                    #print("Raising from {} to {}, cp was {}".format(curr_best_score, tmp_score, cp))
                    break
                    # cp=np.array(cp)*2
                elif tmp_score < curr_best_score:
                    #print("FALLING from {} to {}, cp was {}".format(curr_best_score, tmp_score, cp))
                    curr_best_score = tmp_score
                    projection = tmp_projection
                    change_pro = cp
                    break
                    # cp=np.array(cp)*2
                else:
                    #print ("Score stays at {} in direction {}".format(tmp_score, cp))
                    # Function value did not change. Our stepwidth was to small.
                    cp = np.array(cp) * 2
            else:
                pass  # print ("Projection: No change in direction ", cp)
        if not np.all(change_pro == np.array((0., 0.))):
            # If we didn't change offset, projection or rotation, found still stays True
            found = False
        curr_best_pro = curr_best_pro + change_pro
        if found or not curr_best_score:
            break
    else:
        pass
        #print("Max-iter reached:", maxiter)

    # This is just to detect bugs:
    tmp_projection = fhp.Projection2D(cg, from_polar([1] + list(curr_best_pro)),
                                      project_virtual_atoms=virtual_atoms)
    box = get_box(tmp_projection, scale, curr_best_offs)
    img, _ = tmp_projection.rasterize(dpi, bounding_square=box, warn=False,
                                      rotate=curr_best_rotation)
    best_score = distance(ref_img, img)
    assert best_score == curr_best_score
    box = get_box(projection, scale, curr_best_offs)
    img, _ = projection.rasterize(dpi, bounding_square=box, warn=False,
                                  rotate=curr_best_rotation)
    best_score = distance(ref_img, img)
    assert best_score == curr_best_score
    ######

    return curr_best_score, img, [curr_best_pro, curr_best_rotation, curr_best_offs]


def try_parameters(ref_img, scale, cg,
                   rotations=(0, 180), offsets=(np.array([0, 0]),), proj_directions=None,
                   virtual_atoms=True):
    """
    Try all combinations of the given starting parameters 
    (offset, in-plane rotation and projection direction)
    and find the ones with the shorthest Huasdorff distance.

    :param ref_img: The reference image. A boolean square matrix.
    :param scale: The edge length in Angstrom of the reference image.
    :param cg: The coarse grain RNA to match to the projection.

    :param rotations: A list of in-plane rotations in degrees.
    :param offsets: A list/array of np.arrays of the type np.array([x,y]).
    :param proj_directions: A list of projection_directions (tuples theta, phi) or None.
                            If None, uses cg.project_from.
    :param virtual_atoms: Boolean. If False, do not project virtual atoms (faster)

    :returns: A triple: (best_distance, best_image, best_parameters)
              Where best_distance is a float, best_image a matrix and best params is a triple:
              np.array([theta, phi]), degrees, np.array([x_offset, y_offset])
    """
    if proj_directions is None:
        proj_directions = [to_polar(cg.project_from)[1:]]
    best_distance = float("inf")
    params = [np.array([0, 0]), 0, np.array([0, 0])]
    best_img = None
    for rot in rotations:
        for offset in offsets:
            for proj in proj_directions:
                dist, img, params = locally_minimal_distance(ref_img, scale, cg,  rot, offset, proj,
                                                             maxiter=6, virtual_atoms=virtual_atoms)
                if dist < best_distance:
                    best_distance = dist
                    best_params = params
                    best_img = img
    return best_distance, best_img, best_params


def get_start_points(numPoints):
    """
    Return numPoints equally-distributed points on the unit sphere in polar coordinates.

    Implements the 2nd algorithm from https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
    """
    a = 4 * math.pi / numPoints
    d = math.sqrt(a)
    Mt = int(math.pi / d)
    dt = math.pi / Mt
    df = a / dt
    points = []
    for m in range(Mt):
        theta = math.pi * (m + 0.5) / Mt
        Mf = int(2 * math.pi * math.sin(theta) / df)
        for n in range(Mf):
            phi = 2 * math.pi * n / Mf
            points.append(np.array((theta, phi)))
    return points


def arclength(theta1, phi1, theta2, phi2):
    return (math.acos(math.cos(theta1) * math.cos(theta2) +
                      math.sin(theta1) * math.sin(theta2) * math.cos(phi1 - phi2)))


def get_start_points_near(numPoints, ref_theta, ref_phi):
    a = 4 * math.pi / numPoints
    d = math.sqrt(a)
    Mt = int(math.pi / d)
    dt = math.pi / Mt
    df = a / dt
    points = []
    for m in range(Mt):
        theta = math.pi * (m + 0.5) / Mt
        Mf = int(2 * math.pi * math.sin(theta) / df)
        for n in range(Mf):
            phi = 2 * math.pi * n / Mf
            points.append(np.array((theta, phi)))
    return [p for p in points if arclength(p[0], p[1], ref_theta, ref_phi) < math.pi / 6]


def get_longest_img_diameter(img, scale):
    maxl = 0
    for x1, y1 in np.transpose(np.where(img)):
        for x2, y2 in np.transpose(np.where(img)):
            l = (x1 - x2 + 1)**2 + (y1 - y2 + 1)**2
            if l > maxl:
                maxl = l
    return math.sqrt(maxl) * scale / len(img)


def _try_startpoints(ref_img, scale, cg, start_points, starting_rotations,
                     starting_offsets, local_maxiter, virtual_atoms, use_heuristic, distance, verbose):
    # Longest extention in the image
    longest_distance_image = get_longest_img_diameter(ref_img, scale)
    if longest_distance_image == 0:
        raise ValueError("Reference image probably empty")

    # We count, how often certain heuristics kicked in
    la_heur = 0
    no_heur = 0
    score_heur = 0
    #
    dpi = len(ref_img)
    best_score = float('inf')
    decrease = float("inf")
    for i, project_dir in enumerate(start_points):
        sys.stdout.write("{:2.0%}\r".format(i / len(start_points)))  # Progress
        sys.stdout.flush()
        proj = fhp.Projection2D(cg, from_polar([1] + list(project_dir)),
                                project_virtual_atoms=virtual_atoms)

        if use_heuristic:
            # 4 pixels is arbitrary heuristic
            if abs(proj.longest_axis - longest_distance_image) > 6 * scale / dpi:
                #print("abs({}-{})={}>{}".format(proj.longest_axis, longest_distance_image, abs(proj.longest_axis-longest_distance_image), 6*scale/dpi))
                # proj.plot(show=True)
                la_heur += 1
                continue
        loc_best_rot = 0
        loc_best_offs = np.array([0, 0])
        loc_best_score = float("inf")
        for rot, offset in it.product(starting_rotations, starting_offsets):
            box = get_box(proj, scale, offset)
            img, _ = proj.rasterize(
                dpi, bounding_square=box, warn=False, rotate=rot)
            score = distance(ref_img, img)
            if score < loc_best_score:
                loc_best_score = score
                loc_best_rot = rot
                loc_best_offs = offset
                # loc_best_img=img
        #print("HEUR?", loc_best_score, best_score, best_score+(decrease*1.5))
        if use_heuristic and loc_best_score > best_score + (decrease * 1.75):
            score_heur += 1
            # print("Score")
            # proj.plot(show=True)
            continue
        #import matplotlib.pyplot as plt
        #fig, ax=plt.subplots(3)
        #ax[0].imshow(ref_img, interpolation="none", cmap='gray')
        #ax[1].imshow(loc_best_img, interpolation="none", cmap='gray')
        # ax[0].set_title("Reference")
        #ax[1].set_title("Initial distance {}".format(loc_best_score))
        no_heur += 1
        score, img, params = locally_minimal_distance(ref_img, scale, cg,
                                                      loc_best_rot, loc_best_offs, project_dir,
                                                      maxiter=local_maxiter,
                                                      virtual_atoms=virtual_atoms, distance=distance)
        #ax[2].imshow(img, interpolation="none", cmap='gray')
        #ax[2].set_title("Optimized to {}".format(score))
        # plt.show()
        if score < best_score:
            best_score = score
            best_img = img
            best_params = params
            decrease = loc_best_score - score  # Can increase or decrease
            #print("Decrease {} to {}".format(loc_best_score, score), decrease)
        if best_score == 0:
            break
    sys.stdout.write("   \r")
    if verbose:
        print("Global optimization performe!\n"
              "\t{} local optimizations performed\n"
              "\t{} skipped (longest distance)\n"
              "\t{} skipped (score)".format(no_heur, la_heur, score_heur))
    if no_heur == 0:
        if verbose:
            print("During global optimization, NO SMALL DISTANCE could be found.\n"
                  "\t{} attempts were skipped because the diameters didn't match.\n".format(la_heur))
        return float("inf"), None, [np.array([0, 0]), 0, np.array([0, 0])]
    return best_score, best_img, best_params


def globally_minimal_distance(ref_img, scale, cg,
                              start_points=40,
                              starting_rotations=(0, 180),
                              starting_offsets=(np.array([0, 0]), ),
                              local_maxiter=5, use_heuristic=True, virtual_atoms=True,
                              verbose=False, distance=hausdorff_distance):
    """
    Global minimization of Hausdorff distance.

    Uses several Heuristics to speed up the process.

    :param ref_img: The reference image. A boolean square matrix.
    :param scale: The edge length in Angstrom of the reference image.
    :param cg: The coarse grain RNA to match to the projection.

    :param start_points: Number of starting projection directions. For each, local
                         optimization with maxiter steps is performed.
    :param starting_rotations: A list of in-plane rotations in degrees.
    :param starting_offsets: A list/array of np.arrays of the type np.array([x,y]).

    :param local_maxiter: Maximal iteration for each local optimization
    :param use_heuristic: A Boolean
    :param virtual_atoms: Boolean. If False, do not project virtual atoms (faster)
    :param verbose: If True, print a summary at the end.

    :returns: A triple: (best_distance, best_image, best_parameters)
              Where best_distance is a float, best_image a matrix and best params is a triple:
              np.array([theta, phi]), degrees, np.array([x_offset, y_offset])
    """

    best_score = float("inf")
    decrease = float("inf")
    sp = get_start_points(start_points)

    # We want to cover many different directions with the first few iterations.
    r = random.Random()
    r.seed(1)
    r.shuffle(sp)
    # random.shuffle(sp)

    best_score, best_img, best_params = _try_startpoints(ref_img, scale, cg, sp,
                                                         starting_rotations, starting_offsets, local_maxiter, virtual_atoms,
                                                         use_heuristic, distance, verbose)
    #import matplotlib.pyplot as plt
    #fig, ax=plt.subplots(2)
    #ax[0].imshow(ref_img, interpolation="none", cmap='gray')
    #ax[1].imshow(best_img, interpolation="none", cmap='gray')
    # ax[0].set_title("Reference")
    #ax[1].set_title("distance {}".format(best_score))
    # plt.show()

    # Global search in vicinity of best match
    sp = get_start_points_near(
        10 * start_points, best_params[0][0], best_params[0][1])
    print("Current best score {}, refining on {} points".format(best_score, len(sp)))
    best_score1, best_img1, best_params1 = _try_startpoints(ref_img, scale, cg, sp,
                                                            [best_params[1]], [
                                                                best_params[2]], local_maxiter, virtual_atoms,
                                                            use_heuristic, distance, verbose)
    if best_score1 < best_score:
        best_score = best_score1
        best_params = best_params1
        best_img = best_img1
        print("Refinement helped.", best_score)
    #import matplotlib.pyplot as plt
    #fig, ax=plt.subplots(2)
    #ax[0].imshow(ref_img, interpolation="none", cmap='gray')
    #ax[1].imshow(best_img, interpolation="none", cmap='gray')
    # ax[0].set_title("Reference")
    #ax[1].set_title("distance {}".format(best_score))
    # plt.show()

    # Refinement of the best match
    score, img, params = locally_minimal_distance(ref_img, scale, cg,
                                                  best_params[1], best_params[2], best_params[0],
                                                  maxiter=20 * local_maxiter,
                                                  advanced=True,
                                                  virtual_atoms=virtual_atoms, distance=distance)

    return score, img, params
