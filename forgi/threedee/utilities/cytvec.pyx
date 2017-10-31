# cython: profile=True
import numpy as np
cimport cython
cimport numpy as np
from libc.math cimport sqrt, cos, sin

@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double magnitude(np.ndarray[double, ndim=1] vec) except -1:
    return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])

@cython.wraparound(False)  # turn off negative index wrapping for entire function
def create_orthonormal_basis(np.ndarray[double, ndim=1] vec1, np.ndarray[double, ndim=1] vec2):
    vec1/=magnitude(vec1)
    vec2/=magnitude(vec2)
    cdef np.ndarray[double, ndim=1] vec3 = np.array([vec1[1]*vec2[2]-vec1[2]*vec2[1],  vec1[2]*vec2[0]-vec1[0]*vec2[2], vec1[0]*vec2[1]-vec1[1]*vec2[0] ])
    return np.array([vec1, vec2, vec3])
