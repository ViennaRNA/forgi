# cython: profile=True
import numpy as np
cimport cython
cimport numpy as np
from libc.math cimport sqrt, cos, sin, M_PI, acos


cdef class Vector:
    cdef readonly double x,y,z
    def __init__(self):
        self.x=0.
        self.y=0.
        self.z=0.
    @staticmethod
    cdef Vector from_buffer(double[:] input):
        cdef Vector out = Vector()
        out.x=input[0]
        out.y=input[1]
        out.z=input[2]
        return out

    cdef Vector add(self, Vector other):
        cdef Vector out = Vector()
        out.x = self.x + other.x
        out.y = self.y + other.y
        out.z = self.z + other.z
        return out

    cdef void times(self, double n):
        self.x*=n
        self.y*=n
        self.z*=n

    cdef double square(self):
        return self.x**2+self.y**2+self.z**2

    def __str__(self):
        return "{},{},{}".format(self.coords[0], self.coords[1], self.coords[2])

cdef class Matrix:
    cdef readonly double x0,x1,x2,y0,y1,y2,z0,z1,z2
    def __init__(self):
        pass
    @staticmethod
    cdef Matrix from_buffer(self, double[:,:] input):
        self.x0 = input[0][0]
        self.x1 = input[0][1]
        self.x2 = input[0][2]
        self.y0 = input[1][0]
        self.y1 = input[1][1]
        self.y2 = input[1][2]
        self.z0 = input[2][0]
        self.z1 = input[2][1]
        self.z2 = input[2][2]

def test_vec(a,b):
    cdef Vector avec = Vector.from_buffer(a)
    cdef Vector bvec = Vector.from_buffer(b)
    cdef Vector cvec = avec.add(bvec)
    cdef Vector dvec = cvec.add(avec)
    return dvec

def test_np(a,b):
    c=a+b
    d=c+a
    return d
