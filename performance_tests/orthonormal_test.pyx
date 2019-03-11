# cython: profile=True
import numpy as np
cimport cython
cimport numpy as np
from libc.math cimport sqrt, cos, sin, M_PI, acos

@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double magnitude(double [:] vec ) except -1:
    return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])


@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef np.ndarray[double, ndim=2] create_orthonormal_basis2(
                             np.ndarray[double, ndim=1] vec1,
                             np.ndarray[double, ndim=1] vec2): # Optimized
    cdef double ma, mb, mc
    cdef double ax,ay,az,bx,by,bz
    ma =magnitude(vec1)
    mb =magnitude(vec1)
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
cpdef np.ndarray[double, ndim=2] create_orthonormal_basis3(
                             double[:] vec1,
                             double[:] vec2): # Optimized
    cdef double ma, mb, mc
    cdef double ax,ay,az,bx,by,bz
    ma =magnitude(vec1)
    mb =magnitude(vec1)
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
cpdef np.ndarray[double, ndim=2] create_orthonormal_basis(
                        np.ndarray[double, ndim=1] vec1,
                        np.ndarray[double, ndim=1] vec2): # Slower
    vec1/=magnitude(vec1)
    vec2/=magnitude(vec2)
    cdef np.ndarray[double, ndim=1] vec3 = np.array([vec1[1]*vec2[2]-vec1[2]*vec2[1],  vec1[2]*vec2[0]-vec1[0]*vec2[2], vec1[0]*vec2[1]-vec1[1]*vec2[0] ])
    vec3/=magnitude(vec3)
    return np.array([vec1, vec2, vec3])
