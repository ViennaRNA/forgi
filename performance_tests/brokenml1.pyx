# cython: profile=True
import numpy as np
cimport cython
cimport numpy as np
from libc.math cimport sqrt, cos, sin, M_PI, acos

@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double magnitude(double [:] vec) except -1:
    return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])


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


@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef np.ndarray[double, ndim=1] transform_coords(np.ndarray[double, ndim=2] transposed_stem1_basis,
                           float r, float u, float v):
    cdef np.ndarray[double, ndim=1] xyz = spherical_polar_to_cartesian(r,u,v)
    cdef np.ndarray[double, ndim=1] stem2_start = np.dot(transposed_stem1_basis, xyz)
    return stem2_start

cdef np.ndarray[double, ndim=2] rotation_matrix(char axis, double theta):
    cdef double c = cos(theta)
    cdef double s = sin(theta)
    if axis=="y":
        return np.array([[ c,   0., -s],
                         [ 0.,  1.,  0.],
                         [ s,   0.,  c]])
    elif axis=="z":
        return np.array([[ c,   s,  0.],
                         [-s,   c,  0.],
                         [ 0.,  0., 1.]])
    else:
        assert False

@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef np.ndarray[double, ndim=1] twist2_orient_from_stem1_1(np.ndarray[double, ndim=2] stem1_basis, double u, double v, double t):
    cdef np.ndarray[double, ndim=1] twist2_new = np.array([0., cos(t), sin(t)])

    cdef np.ndarray[double, ndim=2] rot_mat1 = rotation_matrix("z", v)
    cdef np.ndarray[double, ndim=2] rot_mat2 = rotation_matrix("y", u - M_PI / 2.)

    cdef np.ndarray[double, ndim=2] rot_mat = np.dot(rot_mat2, rot_mat1)
    #assert np.allclose(nl.inv(rot_mat), rot_mat.T)
    twist2_new = np.dot(rot_mat.T, twist2_new)

    cdef np.ndarray[double, ndim=1] twist2_new_basis = np.dot(stem1_basis, twist2_new)

    return twist2_new_basis

@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef np.ndarray[double, ndim=1] spherical_polar_to_cartesian(double r, double u, double v):
    cdef double x=r*sin(u)*cos(v)
    cdef double y=r*sin(u)*sin(v)
    cdef double z=r*cos(u)
    return np.array([x,y,z])


def get_fixed_stemvec(cg, broken_ml_name, fixed_stem_name):
    sides = cg.get_sides(fixed_stem_name, broken_ml_name)
    cdef np.ndarray[double, ndim=1] fixed_s_vec = cg.coords.get_direction(fixed_stem_name)
    if sides[0]==0:
        fixed_s_vec = -fixed_s_vec
    cdef np.ndarray[double, ndim=1] fixed_s_pos = cg.coords[fixed_stem_name][sides[0]]
    cdef np.ndarray[double, ndim=1] s_twist = cg.twists[fixed_stem_name][sides[0]]

    return fixed_s_pos, fixed_s_vec, s_twist

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

def get_broken_ml_deviation(cg, broken_ml_name, fixed_stem_name, virtual_stat):
    cdef np.ndarray[double, ndim=1] fixed_s_pos, fixed_s_vec, s_twist
    cdef np.ndarray[double, ndim=1] orig_coords0, orig_coords1, orig_twist

    fixed_s_pos, fixed_s_vec, s_twist = get_fixed_stemvec(cg, broken_ml_name, fixed_stem_name)
    orig_coords0, orig_coords1, orig_twist = get_orig_coords(cg, broken_ml_name, fixed_stem_name)

    cdef np.ndarray[double, ndim=1] v_stat = np.asarray([virtual_stat.r1,
                                                             virtual_stat.u1,
                                                             virtual_stat.v1,
                                                             virtual_stat.u,
                                                             virtual_stat.v,
                                                             virtual_stat.t])
    return _get_broken_ml_dev_core(fixed_s_pos, fixed_s_vec, orig_coords0, orig_coords1, orig_twist, s_twist, v_stat)


@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.wraparound(False)
cdef _get_broken_ml_dev_core(np.ndarray[double, ndim=1] fixed_s_pos,
                             np.ndarray[double, ndim=1] fixed_s_vec,
                             np.ndarray[double, ndim=1] orig_coords0,
                             np.ndarray[double, ndim=1] orig_coords1,
                             np.ndarray[double, ndim=1] orig_twist,
                             np.ndarray[double, ndim=1] s_twist,
                             np.ndarray[double, ndim=1] virtual_stat):
    cdef np.ndarray[double, ndim=2] fixed_stem_basis = create_orthonormal_basis(fixed_s_vec, s_twist)
    cdef np.ndarray[double, ndim=2] transposed_stem1_basis = fixed_stem_basis.transpose()

    cdef np.ndarray[double, ndim=1] vbulge_vec = transform_coords(transposed_stem1_basis,
                                                  virtual_stat[0],
                                                  virtual_stat[1],
                                                  virtual_stat[2])
    cdef np.ndarray[double, ndim=1] vstem_vec = transform_coords(transposed_stem1_basis,
                                                  1.,
                                                  virtual_stat[3],
                                                  virtual_stat[4])
    cdef np.ndarray[double, ndim=1] vstem_twist = twist2_orient_from_stem1_1(transposed_stem1_basis,
                                                   virtual_stat[3],
                                                   virtual_stat[4],
                                                   virtual_stat[5])

    vstem_vec *=5
    cdef np.ndarray[double, ndim=1] vstem_coords0 = fixed_s_pos + vbulge_vec
    cdef np.ndarray[double, ndim=1] vstem_coords1 = vstem_coords0 + vstem_vec

    cdef np.ndarray[double, ndim=1] orig_stem_vec = orig_coords1 - orig_coords0
    cdef np.ndarray[double, ndim=1] true_bulge_vec = orig_coords0 - fixed_s_pos

    cdef double pos_dev = ( vec_distance(orig_coords0, vstem_coords0) )
    cdef double ang_dev = vec_angle(vstem_vec, orig_stem_vec)
    cdef double twist_dev = vec_angle(orig_twist, vstem_twist)
    return pos_dev, ang_dev, twist_dev

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double vec_angle(np.ndarray[double, ndim=1] vec1,
                                    np.ndarray[double, ndim=1] vec2):
        '''
        Get the angle between two vectors using the identity:

        `A * B = |A||B| cos t`

        Where A and B are two vectors and t is the angle between themath.

        :param vec1: The first vector (A)
        :param vec2: The second vector (B)
        :return: The angle between the two vectors.
        '''

        vec1 = normalize(vec1)
        vec2 = normalize(vec2)

        cdef double d = np.dot(vec1, vec2)

        # this shouldn't happen, but sometimes it does, presumably because
        # of rounding errors
        if d >= 1.:
            d = 1.
        if d <= -1.:
            d = -1.

        return  acos(d)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray[double, ndim=1] normalize(np.ndarray[double, ndim=1] vec):
    '''
    Normalize a vector so that its magnitude becomes 1.0 while
    its direction remains the same.

    :param vec: The vector in question.
    :return: A normalized version of the vector.
    '''
    cdef double mag = magnitude(vec)
    return vec / mag

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double vec_distance( np.ndarray[double, ndim=1] vec1,
                                         np.ndarray[double, ndim=1] vec2):
    cdef np.ndarray[double, ndim=1] direction = vec2-vec1
    return sqrt(np.dot(direction, direction))
