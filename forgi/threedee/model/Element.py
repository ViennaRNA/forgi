from __future__ import absolute_import, division, print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.utils import viewkeys
import numpy as np
from collections import Mapping
import logging
log=logging.getLogger(__name__)

class CoordinateStorage(Mapping):
    """
    Provides a dictionary-like interface for the access to coordinates of elements,
    while the coordinates are actually stored as numpy array to allow for manipultion
    with array operations.
    """
    def __init__(self, element_names, on_change = lambda key: None):
        """
        :param defines: The keys that will be used.
        :param on_change: A function that will be called, whenever coordinates are changed.
                        It will receive the key for which the coordinates are changed as sole argument.
        """
        #Initialize coordinate array to NANs
        self._dimensions = 3
        self._coords_per_key = 2 
        self._coordinates = np.ones((self._coords_per_key*len(element_names),self._dimensions))*np.nan
        self._elem_names = { elem: position for position, elem in enumerate(element_names)}
            
        #: on-change function is called whenever coordinates are modified.
        self.on_change = on_change
    def _indices_for(self, elem_name):
        try:
            i = self._elem_names[elem_name]
        except ValueError:
            raise KeyError("Invalid index {}".format(elem_name))
        ret = []
        for j in range(self._coords_per_key):
            ret.append(2*i+j)
        return ret
    def __getitem__(self, elem_name):
        """
        Get the coordinates for a certain coarse grained elements.

        :param elem_name: Either a string with the element name ("s1", "m0", ...)
                          or a sequence of element names.
        :returns: A copy of the corresponding coordinates.
                  If elem_name is a single string, return a tuple of coordinates.
                  If elem_name is a sequence of strings, return a 2*len(elem_name)x3 numpy array.
        """
        try: #Single element name
            indices = self._indices_for(elem_name)
            return tuple(self._coordinates[i] for i in indices)
        except TypeError: # Sequence of element names or invalid
            #Assume sequence of elements
            indices = []
            try:
                for elem in elem_name:
                    try:
                        indices+=self._indices_for(elem)
                    except KeyError as e:
                        raise KeyError("Invalid index: Neither '{}' nor '{}' are valid element "
                                       "names.".format(elem_name, elem_name[0]))
            except TypeError: #elem_name is not iterable
                raise KeyError("Invalid index: Indices must be of type string or sequence of strings. Found {}".format(elem_name))
            return self._coordinates[indices] #Advanced numpy indexing yields a copy.
    def get_direction(self, elem_name): #This assumes the stored coordinates are points not directions 
        assert self._coords_per_key == 2 #Or else a direction does not make sense
        indices = self._indices_for(elem_name)
        return self._coordinates[indices[1]]-self._coordinates[indices[0]]
    
    def __setitem__ (self, key, value):
        if len(value)!=self._coords_per_key:
            raise ValueError("Value must be a {}-tuple of coordinates".format(self._coords_per_key))
        for i in range(self._coords_per_key):
            if len(value[i])!=self._dimensions:
                raise ValueError("Coordinates must have {} dimensions, "
                                 "found {} for entry {}".format(self._dimensions, len(value[i]), i))
        indices = self._indices_for(key)
        for i, index in enumerate(indices):
            self._coordinates[index] = value[i]
        self.on_change(key)
    def __contains__(self, key):
        return key in self._elem_names
    def __iter__(self):
        return iter(self._elem_names)
    def __len__(self):
        return len(self._elem_names)
    def rotate(self, rotation_matrix):
        """
        Rotate all coordinates using the given rotation matrix.
        """
        rotation_matrix = np.asarray(rotation_matrix)
        if rotation_matrix.shape != (3,3):            
            raise ValueError("Rotation matrix does not have the correct shape!")
        self._coordinates = np.dot(self._coordinates, rotation_matrix.T)
        for key in self._elem_names:
            self.on_change(key)
    def get_array(self):
        return np.copy(self._coordinates)

    def __str__(self):
        lines=[]
        for elem in self._elem_names:
            lines.append("{}: ".format(elem)+",  ".join(str(val) for val in self[elem]))
        return("\n".join(lines))
    
    def __eq__(self, other):
        log.debug("Testing equality")
        if type(self)!=type(other):
            log.debug("Typecheck failed")
            return NotImplemented                   
        if viewkeys(self._elem_names) != viewkeys(other._elem_names):
            log.debug("Keys different: self only: {}, other only: {}".format(
                                        viewkeys(self._elem_names)-viewkeys(other._elem_names), 
                                        viewkeys(other._elem_names) - viewkeys(self._elem_names)))
            return False
        if np.all(np.isnan(self._coordinates)) and np.all(np.isnan(other._coordinates)):
            log.debug("True: All is NAN")
            return True
        for key in self:
            if not np.allclose(self[key], other[key]):
                log.debug("Values for key {} different".format(key))
                return False
        return True
    
    def __ne__(self, other):
        return not self == other