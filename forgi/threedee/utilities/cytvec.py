#!/usr/bin/python

#import timeit, sys

import numpy as np
import math as m

def magnitude(vec):
    '''
    Return the magnitude of a vector (|V|).

    @param vec: The vector in question.
    @return: The magnitude of the vector.
    '''
    #return ftuc.magnitude(vec)
    return m.sqrt(np.dot(vec, vec))

def vec_distance(vec1, vec2):
    #return ftuc.vec_distance(vec1, vec2)
    return m.sqrt(np.dot(vec2 - vec1, vec2 - vec1))
