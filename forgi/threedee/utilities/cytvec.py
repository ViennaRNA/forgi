#!/usr/bin/python

#import timeit, sys

import numpy as np
import math as m

def magnitude(vec):
    '''
    Return the magnitude of a vector `(|V|)`.

    :param vec: The vector in question.
    :return: The magnitude of the vector.
    '''
    return m.sqrt(np.dot(vec, vec))

def vec_distance(vec1, vec2):
    if not isinstance(vec1, np.ndarray):
        vec1=np.array(vec1)
        vec2=np.array(vec2)
    return m.sqrt(np.dot(vec2 - vec1, vec2 - vec1))
