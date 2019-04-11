# cython: profile=True
# distutils: language = c++
# distutils: sources = forgi/threedee/utilities/broken_ml_core.cpp

import numpy as np
cimport cython
cimport numpy as np
from libc.math cimport sqrt, cos, sin, acos
from libcpp.vector cimport vector


from forgi.threedee.utilities.broken_ml_core_py cimport Vect, _get_broken_ml_dev_core, testdot

@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double magnitude(double [:] vec) except -1:
    return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])

@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef np.ndarray[double, ndim=2] transposed_inverted(
                             np.ndarray[double, ndim=2] b):
    cdef double det = ( b[0,0]*(b[1,1]*b[2,2]-b[1,2]*b[2,1])
                      - b[0,1]*(b[1,0]*b[2,2]-b[1,2]*b[2,0])
                      + b[0,2]*(b[1,0]*b[2,1]-b[1,1]*b[2,0]))
    return np.array([
                        [(b[1,1]*b[2,2]-b[1,2]*b[2,1])/det, -(b[1,0]*b[2,2]-b[2,0]*b[1,2])/det, (b[1,0]*b[2,1]-b[2,0]*b[1,1])/det],
                        [-(b[0,1]*b[2,2]-b[2,1]*b[0,2])/det, (b[0,0]*b[2,2]-b[2,0]*b[0,2])/det, -(b[0,0]*b[2,1]-b[2,0]*b[0,1])/det],
                        [ (b[0,1]*b[1,2]-b[1,1]*b[0,2])/det, -(b[0,0]*b[1,2]-b[1,0]*b[0,2])/det, (b[0,0]*b[1,1]-b[1,0]*b[0,1])/det]
                    ])

@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef np.ndarray[double, ndim=2] create_orthonormal_basis(
                             double[:] vec1,
                             double[:] vec2): # Optimized
    cdef double ma, mb, mc
    cdef double ax,ay,az,bx,by,bz
    ma =magnitude(vec1)
    mb =magnitude(vec2)
    ax=vec1[0]/ma
    ay=vec1[1]/ma
    az=vec1[2]/ma
    bx=vec2[0]/mb
    by=vec2[1]/mb
    bz=vec2[2]/mb

    cdef np.ndarray[double, ndim=2] out = np.zeros((3,3))
    out[0,0]=ax
    out[0,1]=ay
    out[0,2]=az
    out[1,0]=bx
    out[1,1]=by
    out[1,2]=bz
    out[2,0]=ay*bz-az*by
    out[2,1]=az*bx-ax*bz
    out[2,2]=ax*by-ay*bx
    return out



def get_orig_coords(cg, broken_ml_name, fixed_stem_name):
    s1, s2 = cg.edges[broken_ml_name]
    if s1 == fixed_stem_name:
        orig_stem_name = s2
    elif s2 ==fixed_stem_name:
        orig_stem_name = s1
    else:
        raise ValueError("fixed stem {} is not attached to ml {} with "
                         "edges {}".format(fixed_stem_name, broken_ml_name, [s1,s2]))
    sides2 = cg.get_sides(orig_stem_name, broken_ml_name)
    cdef np.ndarray[double, ndim=1] orig_coords0 = cg.coords[orig_stem_name][sides2[0]]
    cdef np.ndarray[double, ndim=1] orig_coords1 = cg.coords[orig_stem_name][sides2[1]]
    cdef np.ndarray[double, ndim=1] orig_twist   =  cg.twists[orig_stem_name][sides2[0]]
    return orig_coords0, orig_coords1, orig_twist

def get_fixed_stemvec(cg, broken_ml_name, fixed_stem_name):
    sides = cg.get_sides(fixed_stem_name, broken_ml_name)
    cdef np.ndarray[double, ndim=1] fixed_s_vec = cg.coords.get_direction(fixed_stem_name)
    if sides[0]==0:
        fixed_s_vec = -fixed_s_vec
    cdef np.ndarray[double, ndim=1] fixed_s_pos = cg.coords[fixed_stem_name][sides[0]]
    cdef np.ndarray[double, ndim=1] s_twist = cg.twists[fixed_stem_name][sides[0]]

    return fixed_s_pos, fixed_s_vec, s_twist

def get_broken_ml_deviation(cg, broken_ml_name, fixed_stem_name, virtual_stat):
    cdef double[:] fixed_s_pos, fixed_s_vec, s_twist
    cdef double[:] orig_coords0, orig_coords1, orig_twist

    fixed_s_pos, fixed_s_vec, s_twist = get_fixed_stemvec(cg, broken_ml_name, fixed_stem_name)
    orig_coords0, orig_coords1, orig_twist = get_orig_coords(cg, broken_ml_name, fixed_stem_name)

    cdef Vect a =  _get_broken_ml_dev_core(
                Vect(fixed_s_pos[0], fixed_s_pos[1], fixed_s_pos[2]),
                Vect(fixed_s_vec[0], fixed_s_vec[1], fixed_s_vec[2]),
                Vect(orig_coords0[0], orig_coords0[1], orig_coords0[2]),
                Vect(orig_coords1[0], orig_coords1[1], orig_coords1[2]),
                Vect(orig_twist[0], orig_twist[1], orig_twist[2]),
                Vect(s_twist[0], s_twist[1], s_twist[2]),
                virtual_stat.r1, virtual_stat.u1, virtual_stat.v1,
                virtual_stat.u, virtual_stat.v, virtual_stat.t)
    return a.x, a.y, a.z
