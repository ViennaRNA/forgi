
cimport numpy as np

ctypedef np.double_t DTYPE_t
#print dir(np)
#from math import sin, cos, sqrt

from libc.math cimport sin,cos, sqrt, acos, atan2, pow

def vec_distance(np.ndarray[DTYPE_t, ndim=1] vec1, np.ndarray[DTYPE_t, ndim=1] vec2):
    cdef double d0 = vec2[0] - vec1[0]
    cdef double d1 = vec2[1] - vec1[1]
    cdef double d2 = vec2[2] - vec1[2]
    
    return sqrt(d0 * d0 + d1*d1 + d2*d2)
    
def magnitude(np.ndarray[DTYPE_t, ndim=1] vec):
    cdef double x = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2])

    return x



