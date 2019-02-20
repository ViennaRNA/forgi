#!/usr/bin/python
from __future__ import division
from builtins import map
from builtins import range
import timeit
import sys
import itertools
from collections import Counter

import forgi.utilities.debug as fud
import numpy as np
import math
import numpy.linalg as nl
import numpy.testing as nt
import random as rand
import warnings
import logging
#import scipy as sp

log = logging.getLogger(__name__)
try:
    profile
except:
    def profile(f):
        return f

# Set this module-level variable to False to disable all asserts in this module
USE_ASSERTS = __debug__

null_array = np.array([0., 0., 0.])

x_array = np.array([1., 0., 0.])
y_array = np.array([0., 1., 0.])
z_array = np.array([0., 0., 1.])

# identity matrix
identity_matrix = np.array([x_array, y_array, z_array])

standard_basis = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
tau = 2 * math.pi

# Seems to be unused


def get_inter_distances(vecs):
    '''
    Calculate all of the distances between the points of vecs.

    :param vecs: a list of vectors (=points)
    :return: a list containing all distances between any two vectors in vecs.
    '''
    distances = []
    for i in range(len(vecs)):
        for j in range(i + 1, len(vecs)):
            distances += [vec_distance(vecs[i], vecs[j])]

    return distances


def get_random_vector(mult=1.):
    """
    Returns a random vector.

    :param mult: Stretch the random vector by this value. This is the longest value allowed for the total length.
    :return: A random vector
    """
    # Using rejection sampling to generate uniform distribution from points inside a sphere.
    # Thanks to http://stats.stackexchange.com/a/7984/90399
    while True:
        vec = np.array([mult * rand.uniform(-1, 1), mult *
                        rand.uniform(-1, 1), mult * rand.uniform(-1, 1)])
        if magnitude(vec) <= mult:
            return vec


def get_orthogonal_unit_vector(vec):
    '''
    Return a vector orthogonal to vec.

    .. note::

        To create a basis, use create_orthonormal_basis instead!
    '''
    vec2 = get_non_colinear_unit_vector(vec)
    vec3 = np.cross(vec, vec2)
    return normalize(vec3)


# def get_random_vector_pair(angle=rand.uniform(0, math.pi)) -> Removed, because it was never used.


def get_double_alignment_matrix(vp1, vp2):
    '''
    Align two sets of two vectors onto each other.

    :param vp1: A pair of two vectors.
    :param vp2: Another pair of two vectors.
    '''
    angle1 = vec_angle(vp1[0], vp1[1])
    angle2 = vec_angle(vp2[0], vp2[1])

    if USE_ASSERTS:
        nt.assert_allclose(angle1, angle2, rtol=1e-7, atol=1e-7)

    # Align the first two segments
    mat1 = get_alignment_matrix(vp1[0], vp2[0])

    # See where the second segment of the second set ends up
    # after the first alignment
    new_vp2_1 = np.dot(mat1, vp2[1])

    comp1 = np.cross(vp1[1], vp1[0])
    #comp1 = np.cross(vp1[0], vp1[1])
    # should be along the plane of vp1[0] and vp1[1]
    comp2 = np.cross(vp1[0], comp1)

    basis1 = create_orthonormal_basis(normalize(vp1[0]), normalize(comp2))
    rej2 = change_basis(new_vp2_1, basis1, standard_basis)

    angle = math.atan2(rej2[2], rej2[1])

    mat2 = rotation_matrix(vp1[0], angle)

    # return np.dot(mat1, mat2)
    return np.dot(mat2, mat1)


def get_alignment_matrix(vec1, vec2):
    '''
    Return a rotation matrix that will align vec1 along vec2.

    :param vec1: The target vector.
    :param vec2: The vector to be aligned.
    '''

    comp = np.cross(vec1, vec2)
    angle = vec_angle(vec1, vec2)

    return rotation_matrix(comp, angle)


def get_non_colinear_unit_vector(vec):
    '''
    Get a unit vector that does not lie on the line defined by vec.

    This is done by creating a vector along the least represented axis in vec.

    :param vec: The vector under consideration.
    :return: A vector along an axis.
    '''
    absvec = [abs(v) for v in vec]
    m = min(absvec)
    ind = absvec.index(m)
    unit = [0., 0., 0.]
    unit[ind] = 1.

    return np.array(unit)


def is_almost_parallel(vec1, vec2):
    """
    Are vec1 and vec2 almost parallel?
    Nothing is parallel to the zero vector!

    :returns: 1 if vec1 and vec2 are parallel,
             -1 if they are antiparallel,
             0 if they are neither.
    """
    CUTOFF = 10**-7
    vec2_clean = []
    for c in vec2:
        if abs(c) < CUTOFF:
            vec2_clean.append(float("nan"))
        else:
            vec2_clean.append(c)
    factors = np.asarray(vec1) / np.asarray(vec2_clean)
    factors = [f for f in factors if not np.isnan(f)]
    log.debug("vec1 %s, vec2 %s, fac %s", vec1, vec2, factors)
    if not factors:
        return 0  # vec2~[0,0,0]
    elif all(np.sign(factors) == np.sign(factors[0])) and all(abs(f - factors[0]) < CUTOFF for f in factors):
        log.debug("Signs: %s", np.sign(factors))
        return np.sign(factors[0])  # returns 0, if vec1==[0,0,0]
    else:
        return 0  # Not (anti-)parallel


'''
# OLD:
def is_almost_parallel(vec1, vec2):
    """
    Returns true, if two vectors are almost parallel

    Note that every vector is parallel to the zero vector.
    """
    CUTOFF=10**-9
    for i in range(len(vec1)):
        if abs(vec2[i])>CUTOFF:
            factor=vec1[i]/vec2[i]
            break
    else:
        log.debug("vec2 is zero-vector")
        return True # vec2 is Zero-vector
    log.debug("Factor is %f", factor)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        log.debug("vec1/vec2 = %s", [abs(vec1[i]/vec2[i]) for i in range(len(vec1))])
    return all((abs(vec1[i])<CUTOFF and abs(vec2[i])<CUTOFF) or (abs(vec2[i])>CUTOFF and abs(vec1[i]/vec2[i]-factor)<CUTOFF) for i in range(len(vec1)))
'''


def line_segments_collinearity(segment1, segment2):
    """
    Quantifies, how collinear (according to some measure) two line segments are.

    :param segment1, segment2: Each a tuple of vectors (start, end)
    """
    dir1 = normalize(segment1[1] - segment1[0])
    dir2 = normalize(segment2[1] - segment2[0])

    # Get the average direction of the two vectors, if oriented correctly.
    s1 = dir1 + dir2
    s2 = dir1 - dir2
    if magnitude(s1) > magnitude(s2):
        sum_vec = s1
    else:
        sum_vec = s2
    # Further more, the line should pass through the center of the
    # 4 points defining the line segments.
    points = np.array(segment1 + segment2)
    center = points.mean(axis=0)
    '''
    # Now plot for verification
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    log.info("Fiorst 2 points: X %s Y %s Z %s ", points[:2,0], points[:2,1], points[:2,2])
    plt.plot(points[:2,0], points[:2,1], points[:2,2],"ro-")
    plt.plot(points[2:,0], points[2:,1], points[2:,2],"bo-")

    line = np.array([center-3*sum_vec, center, center+3*sum_vec])
    plt.plot(line[:,0],line[:,1],line[:,2], "g-o")
    ax.set_xlim(-5,5)
    ax.set_ylim(-5,5)
    ax.set_zlim(-5,5)

    plt.show()'''

    # Now calculate the distances of the 4 points from the fitted line
    d_point_line = point_line_distance(points, center, sum_vec)
    assert len(d_point_line) == 4
    # Return an R**2 value, like or linear regression (1=collinear, <1 worse)
    return 1 - sum(d_point_line**2) / ((points - center)**2).sum()

    # Different version via SVD
    centered_points = points - center
    uu, dd, vv = np.linalg.svd(centered_points)  # , full_matrices=False)
    # Now vv[0] is the direction of the target line.
    # Calculate the deviation of the points from the target line.
    distances_from_line = np.linalg.norm(
        np.cross(vv[0], centered_points), axis=1) / magnitude(vv[0])
    r2 = 1 - (np.sum(distances_from_line**2) / np.sum(centered_points**2))
    return r2


def create_orthonormal_basis(vec1, vec2=None, vec3=None):
    '''
    Create an orthonormal basis using the provided vectors.

    If more than one is provided, it must be orthogonal to
    the others.

    '''
    if vec1 is not None and vec2 is not None and vec3 is None and not USE_ASSERTS:
        try:
            from . import cytvec
        except ImportError as e:
            warnings.warn(
                "Extention modules (cython code) not installed, using slower python version")
        else:
            return cytvec.create_orthonormal_basis(vec1, vec2)

    if vec2 is None:
        vec2 = get_non_colinear_unit_vector(vec1)
        vec2 = np.cross(vec1, vec2)
    else:
        if USE_ASSERTS:
            if round(vec_angle(vec2, vec1), 9) != round(math.pi / 2, 9):
                assert False, ("vec2 {} is not normal to vec1 {}! Angle is {} rad ({} degrees)".format(
                    vec2, vec1, vec_angle(vec2, vec1), math.degrees(vec_angle(vec2, vec1))))

    mag_vec1 = magnitude(vec1)
    if mag_vec1 == 0:
        raise ZeroDivisionError("vec 1 {}  has magnitude 0.".format(vec1))
    vec1 = vec1 / mag_vec1

    mag_vec2 = magnitude(vec2)
    if mag_vec2 == 0:
        raise ZeroDivisionError(
            "vec 2 has magnitude 0. vecs are so far {} and {} ".format(vec1, vec2))
    vec2 = vec2 / mag_vec2

    if vec3 is None:
        vec3 = np.cross(vec1, vec2)

    mag_vec3 = magnitude(vec3)
    if mag_vec3 == 0:
        raise ZeroDivisionError("vec 3 has magnitude 0. vecs are {}, {} and {}".format(
            repr(vec1), repr(vec2), repr(vec3)))
    vec3 = vec3 / mag_vec3

    return np.array([vec1, vec2, vec3])


"""
# Code used for comparing the fastes method of creating an orthonormal basis:
def create_orthonormal_basis1(vec1, vec2=None, vec3=None):
    '''
    Create an orthonormal basis using the provided vectors.

    If more than one is provided, it must be orthogonal to
    the others.
    '''
    if vec2 == None:
        vec2 = get_non_colinear_unit_vector(vec1)
    #else:
    #    nt.assert_allclose(np.dot(vec2, vec1), 0., rtol=1e-7, atol=1e-7)


    vec1 = normalize(np.array(vec1))
    vec2 = normalize(np.array(vec2))

    if vec3 == None:
        vec3 = np.cross(vec1, vec2)

    vec3 = normalize(vec3)

    return np.array([vec1, vec2, vec3])

def time_cob1():
    vec1 = get_random_vector()
    vec2 = get_random_vector()

    basis = create_orthonormal_basis1(vec1, vec2)

def time_cob2():
    vec1 = get_random_vector()
    vec2 = get_random_vector()

    basis = create_orthonormal_basis(vec1, vec2)

def time_cob():
    t1 = timeit.Timer("time_cob1()", "from forgi.threedee.utilities.vector import time_cob1")
    t2 = timeit.Timer("time_cob2()", "from forgi.threedee.utilities.vector import time_cob2")

    print "1: ", t1.repeat(number=100000) #prints [3.045403003692627, 3.0388529300689697, 3.0359420776367188]
    print "2: ", t2.repeat(number=100000) #prints [2.7473390102386475, 2.74338698387146, 2.731964111328125]
"""


def spherical_cartesian_to_polar(vec):
    '''
    Return a parameterization of a vector of 3 coordinates:

    x = r sin u cos v
    y = r sin u sin v
    z = r cos u

    0 <= u <= pi
    -pi <= v <= pi

    Where u is the polar angle and v is the azimuth angle.

    :param vec: A vector of 3 cartesian coordinates.
    :param fast: Do not assert correctnes of result.
    :return: (r, u, v)
    '''
    r = magnitude(vec)
    u = math.acos(vec[2] / r)
    v = math.atan2(vec[1], vec[0])

    if USE_ASSERTS:
        nt.assert_allclose(vec[0], r * math.sin(u) *
                           math.cos(v), rtol=1e-7, atol=1e-7)
    return np.array((r, u, v))


def spherical_polar_to_cartesian(vec):
    '''
    Convert spherical polar coordinates to cartesian coordinates:

    See the definition of spherical_cartesian_to_polar.

    :param vec: A vector of the 3 polar coordinates (r, u, v)
    :return: (x, y, z)
    '''
    (r, u, v) = vec

    x = r * math.sin(u) * math.cos(v)
    y = r * math.sin(u) * math.sin(v)
    z = r * math.cos(u)

    return np.array([x, y, z])


def get_standard_basis(dim):
    '''
    Get a standard basis for the given dimension.

    For 2D, this equals [[1.,0.],[0.,1.]]

    :param dim: The dimension of the vector space.
    :return: A vector of vectors that constitute the standard basis.
    '''

    standard_basis = [[0. for j in range(dim)] for i in range(dim)]
    for i in range(dim):
        standard_basis[i][i] = 1.
    standard_basis = np.array(standard_basis)

    return standard_basis


def change_basis(coords, new_basis, old_basis):
    '''
    Change the basis of coordinates to a new basis. In a regular structure
    we have the coordinates in the regular cartesian coordinate system. For helix-helix
    orientations, however, we want to express the coordinates in a coordinate system
    defined by the first helix.

    The new basis will consist of the axis of the first helix, one of its twist elements
    (the one closest to the second vector) and a third vector orthogonal to the previous
    two.

    # http://tutorial.math.lamar.edu/Classes/LinAlg/ChangeOfBasis.aspx

    :param coords: The coordinates to transform (array of n elements).
    :param new_basis: The new basis vectors (n x n matrix)
    :param old_basis: The old basis for the coordinates(n x n matrix)
    :return: The new coordinates according to the new basis
    '''
    #assert(len(coords) == len(new_basis))
    #assert(len(new_basis) == len(old_basis))

    dim = len(coords)
    standard_coords = np.dot(old_basis.transpose(), coords)
    new_coords = nl.solve(new_basis.transpose(), standard_coords)

    return new_coords


def change_basis_vectorized(coords, new_basis, old_basis):
    """
    Change an array of vectors (coords) from old_basis to new_basis.

    :param coords: A array of coordinates to transform.
    """
    standard_coords = np.tensordot(
        old_basis.transpose(), coords, axes=([-1], [1])).T
    standard_to_new = nl.inv(new_basis.transpose())
    new_coords = np.tensordot(
        standard_to_new, standard_coords, axes=([-1], [1])).T
    return new_coords


""" # CODE USED FOR BENCHMARKING
    # change_basis1 is slightly faster or the same (4.61 vs 4.67)
def change_basis1(coords, new_basis, old_basis):
    '''
    '''

    dim = len(coords)
    standard_coords = np.dot(old_basis.transpose(), coords)
    standard_to_new = nl.inv(new_basis.transpose())
    new_coords = np.dot(standard_to_new, standard_coords)

    return new_coords


def change_basis2(coords, new_basis, old_basis):
    '''
    '''

    dim = len(coords)
    standard_coords = np.dot(old_basis.T, coords)
    new_coords = nl.solve(new_basis.T, standard_coords)

    return new_coords

def change_basis1_benchmark():
    coords = get_random_vector(10.)
    basis1 = np.array([get_random_vector(10.) for i in range(3)])
    basis2 = np.array([get_random_vector(10.) for i in range(3)])

    nc = change_basis1(coords, basis1, basis2)

def change_basis2_benchmark():
    coords = get_random_vector(10.)
    basis1 = np.array([get_random_vector(10.) for i in range(3)])
    basis2 = np.array([get_random_vector(10.) for i in range(3)])

    nc = change_basis2(coords, basis1, basis2)

def time_basis1():
    t1 = timeit.Timer("change_basis1_benchmark()","from forgi.threedee.utilities.vector import change_basis1_benchmark")
    print t1.repeat(number=100000)

def time_basis2():
    t2 = timeit.Timer("change_basis2_benchmark()","from forgi.threedee.utilities.vector import change_basis2_benchmark")
    print t2.repeat(number=100000)
"""


def vector_rejection(a, b):
    '''
    Return the vector rejection of a from b. In other words, return the orthogonal
    projection of a onto the plane orthogonal to b.

    :param a: The vector to be projected.
    :param b: The vector defining the normal of the plane.
    :return: The rejection of the vector a from b. (a - (np.dot(a, b) / np.dot(b, b)) * b)
    '''

    n = np.dot(a, b)
    d = np.dot(b, b)
    return a - (n / d) * b


def rotation_matrix(axis, theta):
    #TODO:  Use rotation_matrix_weave code (deleted in the commit that introduced this comment)
    # as a guide how to implement this in cython in the future for speedup.
    '''
    Calculate the rotation matrix for a CLOCKWISE rotation of theta around axis.
    This is in the opposite direction that is usually used.

    Thanks to unutbu on StackOverflow

    http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector

    :param axis: The axis around which to rotate. A np-array with length 3 or
                 one of "x", "y" and "z", where x:=standard_basis[0] etc.
    :param theta: The angle of rotation (in rad)
    :return: A matrix which can be used to perform the given rotation. The coordinates
             need only be multiplied by the matrix. (np.dot(matrix, vec))
    '''
    if isinstance(axis, (np.ndarray, list)):
        axis = normalize(axis)
        b, c, d = -axis * math.sin(theta / 2)
        a = math.cos(theta / 2)
        return np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                         [2 * (b * c + a * d), a * a + c * c -
                          b * b - d * d, 2 * (c * d - a * b)],
                         [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c]])
    else:
        s = math.sin(theta)
        a = math.cos(theta)
        if axis == "y":
            return np.array([[a,  0, -s],
                             [0,  1,  0],
                             [s,  0,  a]])
        elif axis == "z":
            return np.array([[a,  s, 0],
                             [-s,  a, 0],
                             [0,  0, 1]])
        elif axis == "x":
            return np.array([[1,  0,  0],
                             [0,  a, s],
                             [0, -s,  a]])
        else:
            raise TypeError('Axis must be numpy array or one of "x", "y", "z"')


def _vector_set_rmsd(set1, set2):
    '''
    Calculate the not-centered rmsd between two sets of vectors.
    Currently in forgi.threedee.utilities.pdb, but subject to future deprecation.

    :param set1: A matrix
    :param set2: Another matrix.

    :return: The rmsd between the rows of the matrix.
    '''
    rmsd = 0
    count = 0
    for i in range(len(set1)):
        rmsd += magnitude(set2[i] - set1[i]) ** 2
        count += 1
    rmsd /= count
    return math.sqrt(rmsd)


def get_vector_centroid(crds1):
    '''
    Find the centroid of a set of vectors.

    :param crds: A matrix containing all of the vectors.

    :return: The centroid of the rows of the matrix crds.
    '''
    crds1 = np.asarray(crds1)
    centroid = np.sum(crds1, axis=0)
    centroid /= float(len(crds1))

    return centroid
    if crds1.shape[1] == 3:
        centroid1 = np.array([0., 0., 0.])
    else:
        centroid1 = np.array([0., 0.])

    for i in range(len(crds1)):
        centroid1 += crds1[i]

    centroid1 /= float(len(crds1))

    for i in centroid1:
        if math.isnan(i):
            raise ValueError('nan encountered in centroid: {}, len crds1 = {}.'.format(
                centroid1, len(crds1)))

    return centroid1


def center_on_centroid(crds1):
    centroid1 = get_vector_centroid(crds1)

    crds = np.asarray(crds1)
    return crds - centroid1


def magnitude(vec):
    '''
    Return the magnitude of a vector `(|V|)`.

    This is guaranteed for arbitrary dimensions of vec.

    :param vec: The vector in question.
    :return: The magnitude of the vector.
    '''
    # return np.linalg.norm(vec) #A lot of overhead, if used for a single vector
    return np.sqrt(np.dot(vec, vec))


def det3x3(matrix):
    """return the determinant of a 3x3 matrix"""
    positive = matrix[0, 0] * matrix[1, 1] * matrix[2, 2] + matrix[0, 1] * \
        matrix[1, 2] * matrix[2, 0] + \
        matrix[0, 2] * matrix[1, 0] * matrix[2, 1]
    minus = matrix[2, 0] * matrix[1, 1] * matrix[0, 2] + matrix[2, 1] * \
        matrix[1, 2] * matrix[0, 0] + \
        matrix[2, 2] * matrix[1, 0] * matrix[0, 1]
    return positive - minus


"""
def time_mag1():
    vec1 = get_random_vector()

    return math.sqrt(np.dot(vec1, vec1))

def time_mag2():
    vec1 = get_random_vector()

    return math.sqrt(np.dot(vec1, vec1))

def time_mag():
    t1 = timeit.Timer("time_mag1()", "from forgi.utilities.vector import time_mag1")
    t2 = timeit.Timer("time_mag2()", "from forgi.utilities.vector import time_mag2")

    print t1.repeat(number=10000)
    print t2.repeat(number=10000)
"""


def normalize(vec):
    '''
    Normalize a vector so that its magnitude becomes 1.0 while
    its direction remains the same.

    :param vec: The vector in question.
    :return: A normalized version of the vector.
    '''
    mag = magnitude(vec)
    if mag == 0:  # Numpy would return Nan and raise a RuntimeWarning.
        raise ValueError("Cannot normalize zero- vector!")
    return vec / mag


def vec_angle(vec1, vec2):
    '''
    Get the angle between two vectors using the identity:

    `A * B = |A||B| cos t`

    Where A and B are two vectors and t is the angle between themath.

    :param vec1: The first vector (A)
    :param vec2: The second vector (B)
    :return: The angle between the two vectors.
    '''

    vec1n = normalize(vec1)
    vec2n = normalize(vec2)

    d = np.dot(vec1n, vec2n)

    # this shouldn't happen, but sometimes it does, presumably because
    # of rounding errors
    if d >= 1.:
        d = 1.
    if d <= -1.:
        d = -1.

    angle = math.acos(d)
    return angle


def vec_dot(a, b):
    """
    Vector dot product for vectors of length 3.

    For small vectors of length 3 that are represented as lists and not as np.arary,
    this naive python implementation might be faster than the corresponding numpy implementation.

    If a and b are already numpy arrays, the numpy implementation seems to be faster
    (depending an how numpy was compiled)

    """
    # >>> timeit.timeit('ftuv.vec_dot(a,b)', setup='import numpy as np;import forgi.threedee.utilities.vector as ftuv; a=np.array([1,2,3]); b=np.array([7,2,5]);')
    # 0.7863450050354004
    # >>> timeit.timeit('ftuv.vec_dot(a,b)', setup='import numpy as np;import forgi.threedee.utilities.vector as ftuv; a=[1,2,3]; b=[7,2,5];')
    # 0.27366209030151367
    # >>> timeit.timeit('np.dot(a,b)', setup='import numpy as np; a=[1,2,3]; b=[7,2,5]')
    # 1.4038841724395752
    # >>> timeit.timeit('np.dot(a,b)', setup='import numpy as np; a=np.array([1,2,3]); b=np.array([7,2,5])')
    # 0.6194930076599121
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


'''
def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return c
'''


def seg_intersect(line1, line2):
    """
    Intersection of 2 line segments in 2D space (as lists or numpy array-like).
    :param line1: a tuple/list (a1, a2): The first line segment, from a1 to a2
    :param line2: a tuple/list (b1, b2):The 2nd line segment, from b1 to b2
    """
    a1, a2 = line1
    b1, b2 = line2
    a1 = np.array(a1)
    a2 = np.array(a2)
    b1 = np.array(b1)
    b2 = np.array(b2)
    if max(map(len, [a1, a2, b1, b2])) != 2:
        raise ValueError(
            "Expecting only 2-dimensional vectors. Found higher-dimensional vector: a1={}, a2={}, b1={}, b2={}".format(a1, a2, b1, b2))
    if min(map(len, [a1, a2, b1, b2])) != 2:
        raise ValueError(
            "Expecting only 2-dimensional vectors. Found lower-dimensional vector.")
    if (a1 == a2).all() or (b1 == b2).all():
        raise ValueError(
            "Start and end of a line must not be equal! a1={}, a2={}, b1={}, b2={}".format(a1, a2, b1, b2))
    dxa = a2[0] - a1[0]
    dya = a2[1] - a1[1]
    dxb = b2[0] - b1[0]
    dyb = b2[1] - b1[1]
    num = a1[0] * dya - a1[1] * dxa - b1[0] * dya + b1[1] * dxa
    denom = float(dxb * dya - dyb * dxa)

    if denom == 0:
        # parallel or on same lines
        if dxa == 0:
            t1 = (b1[1] - a1[1]) / dya
            t2 = (b2[1] - a1[1]) / dya
            t1test = t1
        else:
            t1 = (b1[0] - a1[0]) / dxa
            t2 = (b2[0] - a1[0]) / dxa
        if dya == 0:
            t1test = t1
        else:
            t1test = (b1[1] - a1[1]) / dya
        if t1 != t1test:
            return []
        # On same line
        if dxa == 0:
            s1 = (a1[1] - b1[1]) / dyb
            s2 = (a2[1] - b1[1]) / dyb
        else:
            s1 = (a1[0] - b1[0]) / dxb
            s2 = (a2[0] - b1[0]) / dxb
        if all(x < 0 or x > 1 for x in [s1, s2, t1, t2]):
            return []
        ts = min(t1, t2)
        te = max(t1, t2)
        toret = []
        if ts < 0:
            toret.append(a1)
        else:
            toret.append(b1)
        if te < 1:
            toret.append(b2)
        else:
            toret.append(a2)
        return toret
    else:
        s = num / denom
        if s >= 0 and s <= 1:
            c = np.array(b1) + s * (np.array(b2) - np.array(b1))
            if dxa != 0:
                t = (c[0] - a1[0]) / dxa
            else:
                t = (c[1] - a1[1]) / dya
            if t >= 0 and t <= 1:
                return [c]
        return []


def point_line_distance(point, line_start, line_dir):
    """
    Calculate the distance between the point and the line
    through line_point with direction line_dir.

    :param point: A point(an array with shape (3,)) or multiple points (shape n,3)
    """
    if np.shape(point) == (3,):
        return magnitude(np.cross(line_dir, (line_start - point))) / magnitude(line_dir)
    else:
        # More than one point
        return np.linalg.norm(np.cross(line_dir, line_start - point), axis=1) / magnitude(line_dir)


def vec_distance(vec1, vec2):
    """
    The (euclidean) distance between two points vec1 and vec2.

    This is guaranteed to work for arbitrary but equal dimensions of vec1 and vec2.

    :param vec1, vec2: A list or np.array of floats
    :returns: A float
    """
    vec1 = np.asarray(vec1)
    vec2 = np.asarray(vec2)
    direction = vec2 - vec1
    return math.sqrt(np.dot(direction, direction))


@profile
def elements_closer_than(s1_p0, s1_p1, s2_p0, s2_p1, distance):
    '''
    Code copied from line_segment_distance, but with optimizations for fast comparison to distance.

    Code shamelessly translated from:
    http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment

    :param s1_p0: The start of the first segment
    :param s1_p1: The end of the first segment

    :param s2_p0: The start of the second segment
    :param s2_p1: The end of the second segment

    :return: True or False
    '''
    # pylint: disable=E1130  # See https://github.com/PyCQA/pylint/issues/2721
    u = s1_p1 - s1_p0
    v = s2_p1 - s2_p0
    w = s1_p0 - s2_p0
    lenw = magnitude(w)
    a = np.dot(u, u)        # always >= 0
    c = np.dot(v, v)        # always >= 0

    if lenw < distance:
        return True
    if lenw > math.sqrt(a) + math.sqrt(c) + distance:
        return False

    b = np.dot(u, v)

    d = np.dot(u, w)
    e = np.dot(v, w)

    D = a * c - b * b       # always >= 0
    sD = D      # sc = sN / sD, default sD = D >= 0
    tD = D      # tc = tN / tD, default tD = D >= 0

    SMALL_NUM = 0.000001

    # compute the line parameters of the two closest points
    if (D < SMALL_NUM):  # the lines are almost parallel
        sN = 0.0        # force using point P0 on segment S1
        sD = 1.0        # to prevent possible division by 0.0 later
        tN = e
        tD = c
    else:                # get the closest points on the infinite lines
        sN = (b * e - c * d)
        tN = (a * e - b * d)
        if (sN < 0.0):      # sc < 0 => the s=0 edge is visible
            sN = 0.0
            tN = e
            tD = c
        elif (sN > sD):  # sc > 1 => the s=1 edge is visible
            sN = sD
            tN = e + b
            tD = c

    if (tN < 0.0):           # tc < 0 => the t=0 edge is visible
        tN = 0.0
        # recompute sc for this edge
        if (-d < 0.0):
            sN = 0.0
        elif (-d > a):
            sN = sD
        else:
            sN = -d
            sD = a
    elif (tN > tD):      # tc > 1 => the t=1 edge is visible
        tN = tD
        # recompute sc for this edge
        if ((-d + b) < 0.0):
            sN = 0
        elif ((-d + b) > a):
            sN = sD
        else:
            sN = (-d + b)
            sD = a

    # finally do the division to get sc and tc
    sc = 0.0 if abs(sN) < SMALL_NUM else sN / sD
    tc = 0.0 if abs(tN) < SMALL_NUM else tN / tD

    # get the difference of the two closest points
    # dP = w + (sc * u) - (tc * v)  # = S1(sc) - S2(tc)

    return vec_distance(s1_p0 + sc * u, s2_p0 + tc * v) < distance


def line_segment_distance(s1_p0, s1_p1, s2_p0, s2_p1):
    '''
    Calculate the two points on each of the segments that are closest to
    each other. The first segment is defined as p1->p2 and the second as
    p3->p4.

    Code shamelessly translated from:
    http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment

    :param s1_p0: The start of the first segment
    :param s1_p1: The end of the first segment

    :param s2_p0: The start of the second segment
    :param s2_p1: The end of the second segment

    :return: A tuple of points (i1,i2) containing the point i1 on s1
        closest to the point i2 on segment s2.
    '''
    # pylint: disable=E1130  # See https://github.com/PyCQA/pylint/issues/2721
    u = s1_p1 - s1_p0
    v = s2_p1 - s2_p0
    w = s1_p0 - s2_p0

    a = np.dot(u, u)        # always >= 0
    b = np.dot(u, v)
    c = np.dot(v, v)        # always >= 0
    d = np.dot(u, w)
    e = np.dot(v, w)

    D = a * c - b * b       # always >= 0
    sD = D      # sc = sN / sD, default sD = D >= 0
    tD = D      # tc = tN / tD, default tD = D >= 0

    SMALL_NUM = 0.000001

    # compute the line parameters of the two closest points
    if (D < SMALL_NUM):  # the lines are almost parallel
        sN = 0.0        # force using point P0 on segment S1
        sD = 1.0        # to prevent possible division by 0.0 later
        tN = e
        tD = c
    else:                # get the closest points on the infinite lines
        sN = (b * e - c * d)
        tN = (a * e - b * d)
        if (sN < 0.0):      # sc < 0 => the s=0 edge is visible
            sN = 0.0
            tN = e
            tD = c
        elif (sN > sD):  # sc > 1 => the s=1 edge is visible
            sN = sD
            tN = e + b
            tD = c

    if (tN < 0.0):           # tc < 0 => the t=0 edge is visible
        tN = 0.0
        # recompute sc for this edge
        if (-d < 0.0):
            sN = 0.0
        elif (-d > a):
            sN = sD
        else:
            sN = -d
            sD = a
    elif (tN > tD):      # tc > 1 => the t=1 edge is visible
        tN = tD
        # recompute sc for this edge
        if ((-d + b) < 0.0):
            sN = 0
        elif ((-d + b) > a):
            sN = sD
        else:
            sN = (-d + b)
            sD = a

    # finally do the division to get sc and tc
    sc = 0.0 if abs(sN) < SMALL_NUM else sN / sD
    tc = 0.0 if abs(tN) < SMALL_NUM else tN / tD

    # get the difference of the two closest points
    # dP = w + (sc * u) - (tc * v)  # = S1(sc) - S2(tc)

    return (s1_p0 + sc * u, s2_p0 + tc * v)


def closest_point_on_seg(seg_a, seg_b, circ_pos):
    '''
    Closest point between a line segment and a point.

    Lifted from:

    http://doswa.com/2009/07/13/circle-segment-intersectioncollision.html
    '''
    if not isinstance(seg_a, np.ndarray):
        seg_a = np.array(seg_a)
        seg_b = np.array(seg_b)
        circ_pos = np.array(circ_pos)
    seg_v = seg_b - seg_a
    pt_v = circ_pos - seg_a
    mag = math.sqrt(sum(seg_v * seg_v))

    if mag <= 0:
        raise ValueError("Invalid segment length")
    seg_v_unit = seg_v / mag
    proj = pt_v.dot(seg_v_unit)
    if proj <= 0:
        return seg_a.copy()
    if proj >= mag:
        return seg_b.copy()
    proj_v = seg_v_unit * proj
    closest = proj_v + seg_a
    return closest


def cylinder_line_intersection(cyl, line, r):
    '''
    Get the points of intersection between a line and a cylinder.

    If they do not intersect, return an empty list. If the line
    touches the cylinder, then return a 2 point list with two identical points.
    If the line crosses the cylinder, then return a list of 2 points.
    '''

    cyl = np.array(cyl)
    line = np.array(line)
    cyl_vec = cyl[1] - cyl[0]
    #line_vec = line[1] - line[0]

    cyl_basis = create_orthonormal_basis(cyl_vec)

    line_t = change_basis(line.T, cyl_basis, standard_basis).T
    cyl_t = change_basis(cyl.T, cyl_basis, standard_basis).T

    #cyl_vec_t = cyl_t[1] - cyl_t[0]
    line_vec_t = line_t[1] - line_t[0]

    line_vec_t_normed = normalize(line_vec_t)

    # get the point on the line closest to the
    # center of the cylinder
    p = closest_point_on_seg(line_t[0][1:], line_t[1][1:],
                             np.array(cyl_t[0][1:]))

    if line_vec_t[1] == 0:
        return []

    # Figure out the x position by determining how far along
    # the y-coordinate of the segment the closest point is
    x = (line_t[0][0] +
         (line_vec_t[0] *
          (p[0] - line_t[0][1]) / line_vec_t[1]))
    v = p - cyl_t[0][1:]
    o = math.sqrt(sum(v * v))
    p = [x, p[0], p[1]]

    if o > r:
        # no intersection
        return []

    t = math.sqrt(r ** 2 - o ** 2)

    i1 = p - t * line_vec_t_normed
    i2 = p + t * line_vec_t_normed

    endpoints_t = [i1, i2]

    endpoints_t = sorted(endpoints_t, key=lambda x: x[0])
    line_t = sorted(line_t, key=lambda x: x[0])
    cyl_t = sorted(cyl_t, key=lambda x: x[0])

    # the start point will be the highest of the top of the cylinder
    # the end of the line, and the point where the line intersects
    # the surface of the cylinder
    start_points = np.array([cyl_t[0][0], line_t[0][0], endpoints_t[0][0]])

    # the end point will be the lowest of the... " " "
    end_points = np.array([cyl_t[1][0], line_t[1][0], endpoints_t[1][0]])

    start_points = sorted(start_points)
    end_points = sorted(end_points)

    start_param = (start_points[-1] - p[0]) / line_vec_t_normed[0]
    end_param = (end_points[0] - p[0]) / line_vec_t_normed[0]

    endpoints_t[0] = p + start_param * line_vec_t_normed
    endpoints_t[1] = p + end_param * line_vec_t_normed

    '''
    endpoints_t[0] = [start_points[-1],
                      endpoints_t[
    endpoints_t[1][0] = end_points[0]
    real_start = start_points[-1]
    real_end = end_points[0]
    '''

    intersects_t = np.array(endpoints_t)

    if intersects_t[0][0] > intersects_t[1][0]:
        # degenerate case, the line is to the side of the cylinder
        return []

    return change_basis(intersects_t.T, standard_basis, cyl_basis).T

# COVERAGE: Not used in ernwin and forgi at least since forgi v0.3


def pin_fits_two_cyl(cyl1, cyl2, cyl_width):
    '''
    If we create a cone that starts at one end of cylinder1, does it
    enclose cylinder2?

    This function approximates an answer by projecting a circle on
    the plane normal to the axis of cylinder1.

    :param cyl1: The coordinates of cylinder1
    :param cyl2: The coordinates of cylinder2
    :param cyl_width: The widths of the two cylinders
    :return: True or False
    '''
    basis = create_orthonormal_basis(cyl1[1] - cyl1[0])
    cyl2 = np.array([cyl2[0] - cyl1[1], cyl2[1] - cyl1[1]])

    cyl2_t = change_basis(cyl2.T, basis, standard_basis).T
    cone_width = cyl_width
    cyl1_len = magnitude(cyl1[1] - cyl1[0])

    # the cone expands
    cone_width_cyl2_start = cyl_width * (cyl1_len + cyl2_t[0][0]) / cyl1_len
    cone_width_cyl2_end = cyl_width * (cyl1_len + cyl2_t[1][0]) / cyl1_len

    cyl2_start_offset = math.sqrt(cyl2_t[0][1] ** 2 + cyl2_t[0][2] ** 2)
    cyl2_end_offset = math.sqrt(cyl2_t[1][1] ** 2 + cyl2_t[1][2] ** 2)

    if cyl2_start_offset > cone_width_cyl2_start:
        return False
    if cyl2_end_offset > cone_width_cyl2_end:
        return False

    return True

# COVERAGE: Not used in ernwin and forgi at least since forgi v0.3


def point_in_cylinder(pt1, pt2, r, testpt):
    '''
    Determine if testpt is within a cylinder of radius r.

    Translated from:

        http://www.flipcode.com/archives/Fast_Point-In-Cylinder_Test.shtml

    :param pt1: The start of the cylinder axis.
    :param pt2: The end of the cylinder axis.
    :param r: The radius of the cylinder.
    :param testpt: The point we are testing.
    :return: True if the point is within the cylinder, False otherwise.
    '''
    d = pt2 - pt1
    lengthsq = np.dot(d, d)
    radius_sq = r * r

    pd = testpt - pt1
    dot = np.dot(pd, d)

    if dot < 0. or dot > lengthsq:
        # beyond the edges of the cylinder
        return False
    else:
        dsq = np.dot(pd, pd) - (dot * dot) / lengthsq

        if (dsq > radius_sq):
            return False
        else:
            return True


def GetPointsEquiAngularlyDistancedOnSphere(numberOfPoints=45):
    """ each point you get will be of form 'x, y, z'; in cartesian coordinates
        eg. the 'l2 distance' from the origion [0., 0., 0.] for each point will be 1.0
        ------------
        converted from:  http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere )
    """
    dlong = math.pi * (3.0 - math.sqrt(5.0))  # ~2.39996323
    dz = 2.0 / numberOfPoints
    long = 0.0
    z = 1.0 - dz / 2.0
    ptsOnSphere = []
    for k in range(0, numberOfPoints):
        r = math.sqrt(1.0 - z * z)
        ptNew = (math.cos(long) * r, math.sin(long) * r, z)
        ptsOnSphere.append(ptNew)
        z = z - dz
        long = long + dlong

    return ptsOnSphere


def sortAlongLine(start, end, points):
    """
    Sort all points in points along the line from start to end.

    :param start: A point
    :param end: A point
    :param points: A list of points

    :returns: A list containing start, end and all elements of points, sorted by the distance from start
    """
    # print start, end, points
    s_points = points + [start, end]
    s_points.sort(key=lambda x: vec_distance(x, start))
    return s_points


def middlepoint(vec1, vec2):
    """The point in the middle between vec1 and vec2."""
    generator = ((x + vec2[i]) / 2.0 for i, x in enumerate(vec1))
    typ = type(vec1)
    if typ == np.ndarray:
        return np.fromiter(generator, float, len(vec1))
    return typ(generator)


def pair_distance_distribution(points, stepsize=1):
    return pair_distance_distribution_vectorized(points, stepsize)
    # The below is left in here as a not vectorized (slow) reference implementation
    dists = Counter()
    for p1, p2 in itertools.combinations(points, r=2):
        d = vec_distance(p1, p2)
        dists[d // stepsize] += 1
    m = max(dists.keys())
    x = np.arange(0, m + 1)
    return x * stepsize, np.array([dists[xi] for xi in x])

def pair_distance_distribution_vectorized(points, stepsize=1):
    dists = Counter()
    p1s = []
    p2s = []
    for p1, p2 in itertools.combinations(points, r=2):
        p1s.append(p1)
        p2s.append(p2)
    diffs=np.array(p1s)-np.array(p2s)
    lengths=np.sqrt(np.sum(diffs*diffs, axis=1))
    lengths = lengths // stepsize
    for d in lengths:
        dists[d] += 1
    m = max(dists.keys())
    x = np.arange(0, m + 1)
    return x * stepsize, np.array([dists[xi] for xi in x])
