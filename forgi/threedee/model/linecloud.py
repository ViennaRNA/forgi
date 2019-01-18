from __future__ import absolute_import, division, print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
import sys
if sys.version_info < (3,):
    from future.utils import viewkeys
else:
    viewkeys = lambda dic, **kwargs: dic.keys(**kwargs)

import numpy as np
from collections import Mapping
import itertools
import logging
import forgi.threedee.utilities.vector as ftuv

log = logging.getLogger(__name__)

try:
    profile
except:
    def profile(f):
        return f


class CoordinateStorage(Mapping):
    """
    Provides a dictionary-like interface for the access to coordinates of elements,
    while the coordinates are actually stored as numpy array to allow for manipultion
    with array operations.
    """

    def __init__(self, element_names, on_change=lambda key: None):
        """
        :param defines: The keys that will be used.
        :param on_change: A function that will be called, whenever coordinates are changed.
                        It will receive the key for which the coordinates are changed as sole argument.
        """
        # Initialize coordinate array to NANs
        self._dimensions = 3
        self._coords_per_key = 2
        log.debug("Initializing CoordinateStorage with Nans")
        self._coordinates = np.ones(
            (self._coords_per_key * len(element_names), self._dimensions)) * np.nan
        self._elem_names = {elem: position for position,
                            elem in enumerate(element_names)}

        #: on-change function is called whenever coordinates are modified.
        self.on_change = on_change

    @profile
    def _indices_for(self, elem_name):
        try:
            i = self._elem_names[elem_name]
        except ValueError:
            raise KeyError("Invalid index {}".format(elem_name))
        if self._coords_per_key == 2:
            return [2 * i, 2 * i + 1]
        else:
            ret = []
            for j in range(self._coords_per_key):
                ret.append(2 * i + j)
            return ret

    @profile
    def __getitem__(self, elem_name):
        """
        Get the coordinates for a certain coarse grained elements.

        :param elem_name: Either a string with the element name ("s1", "m0", ...)
                          or a sequence of element names.
        :returns: A copy of the corresponding coordinates.
                  If elem_name is a single string, return a tuple of coordinates.
                  If elem_name is a sequence of strings, return a 2*len(elem_name)x3 numpy array.
        """
        try:  # Single element name
            indices = self._indices_for(elem_name)
            return tuple(self._coordinates[i] for i in indices)
        except TypeError:  # Sequence of element names or invalid
            # Assume sequence of elements
            indices = []
            try:
                for elem in elem_name:
                    try:
                        indices += self._indices_for(elem)
                    except KeyError as e:
                        raise KeyError("Invalid index: Neither '{}' nor '{}' are valid element "
                                       "names.".format(elem_name, elem_name[0]))
            except TypeError:  # elem_name is not iterable
                raise KeyError(
                    "Invalid index: Indices must be of type string or sequence of strings. Found {}".format(elem_name))
            # Advanced numpy indexing yields a copy.
            return self._coordinates[indices]

    @profile
    def __setitem__(self, key, value):
        # Commented out for speed-gain
        # if len(value)!=self._coords_per_key:
        #    raise ValueError("Value must be a {}-tuple of coordinates".format(self._coords_per_key))
        # for i in range(self._coords_per_key):
        #    if len(value[i])!=self._dimensions:
        #        raise ValueError("Coordinates must have {} dimensions, "
        #                         "found {} for entry {}".format(self._dimensions, len(value[i]), i))
        indices = self._indices_for(key)
        for i, index in enumerate(indices):
            log.debug("Setting coordinate %s(%d) to %s", key, i, value[i])
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
        if rotation_matrix.shape != (3, 3):
            raise ValueError(
                "Rotation matrix does not have the correct shape!")
        self._coordinates = np.dot(self._coordinates, rotation_matrix.T)
        for key in self._elem_names:
            self.on_change(key)

    def get_array(self):
        """
        Return a copy of the underlying numpy array.
        """
        return np.copy(self._coordinates)

    @property
    def is_filled(self):
        """
        Returns True if all coordinates are set to somethiong different than nan.
        """
        return not np.isnan(self._coordinates).any()

    def __str__(self):
        lines = []
        for elem in self._elem_names:
            lines.append("{}: ".format(elem) + ",  ".join(str(val)
                                                          for val in self[elem]))
        return("\n".join(lines))

    def __repr__(self):
        return "<{} object at {} with {}>".format(type(self).__name__, hex(id(self)), str(self).replace("\n", "; "))

    def __eq__(self, other):
        log.debug("Testing equality")
        if type(self) != type(other):
            log.debug("Typecheck failed")
            return NotImplemented
        if viewkeys(self._elem_names) != viewkeys(other._elem_names):
            log.debug("Keys different: self only: {}, other only: {}".format(
                viewkeys(self._elem_names) - viewkeys(other._elem_names),
                viewkeys(other._elem_names) - viewkeys(self._elem_names)))
            return False
        if np.all(np.isnan(self._coordinates)) and np.all(np.isnan(other._coordinates)):
            log.debug("True: All is NAN")
            return True
        for key in self:
            if not np.allclose(self[key], other[key]):
                log.debug("Values for key {} different: {}!={}".format(
                    key, self[key], other[key]))
                return False
        log.debug("Equal!")
        return True

    def __ne__(self, other):
        return not self == other


class LineSegmentStorage(CoordinateStorage):
    _coords_per_key = 2

    def __init__(self, *args, **kwargs):
        super(LineSegmentStorage, self).__init__(*args, **kwargs)
        self._i_to_elem = {i: elem for elem, i in self._elem_names.items()}
        self.is_centered = False
        self.on_change = self._extended_on_change(self.on_change)

    def _extended_on_change(self, f):
        """
        In addition to the user-supplied function, also reset is_centered.
        """
        def fu(key):
            f(key)
            self.is_centered = False
        return fu

    def center(self):
        self._coordinates = ftuv.center_on_centroid(self._coordinates)
        self.is_centered = True

    # This assumes the stored coordinates are points not directions
    def get_direction(self, elem_name):
        assert self._coords_per_key == 2  # Or else a direction does not make sense
        indices = self._indices_for(elem_name)
        return self._coordinates[indices[1]] - self._coordinates[indices[0]]

    @profile
    def elements_closer_than(self, cutoff, ignore=[]):
        """
        :param ignore: A set of tuples (element name-pairs) to ignore
        """
        # See http://stackoverflow.com/a/18994296/5069869 by Fnord
        # Modified to make use of numpy vectorization.
        i_to_elem = self._i_to_elem

        assert self._coords_per_key == 2
        directions = self._coordinates[1::2] - self._coordinates[::2]
        magnitudes = np.linalg.norm(directions, axis=1, keepdims=True)
        normed_directions = directions / magnitudes
        hits = []
        for i, j in itertools.combinations(range(len(self._elem_names)), 2):
            potential_interaction = tuple(sorted((i_to_elem[i], i_to_elem[j])))
            node1, node2 = potential_interaction
            if (node1, node2) in ignore or (node2, node1) in ignore:
                log.debug("Ignoring nodes %s", potential_interaction)
                continue
            else:
                log.debug("Testing closeness of %s", (node1, node2))
            a0 = self._coordinates[2 * i]
            a1 = self._coordinates[2 * i + 1]
            b0 = self._coordinates[2 * j]
            b1 = self._coordinates[2 * j + 1]
            vec_a0b0 = b0 - a0
            len_a0b0 = ftuv.magnitude(vec_a0b0)
            if len_a0b0 < cutoff:
                log.debug("a0b0 already confirms hit for %s",
                          potential_interaction)
                hits.append(potential_interaction)
                continue
            elif len_a0b0 > cutoff + magnitudes[i] + magnitudes[j]:
                log.debug("a0b0 already rules out a hit for %s",
                          potential_interaction)
                continue  # Cannot be closer than cutoff!
            a_normed = normed_directions[i]
            b_normed = normed_directions[j]
            cross = np.cross(a_normed, b_normed)
            denom = np.sum(cross**2)  # =norm**2
            SMALL_NUMBER = 0.000001
            if denom < SMALL_NUMBER:  # lines are parallel
                d0 = np.dot(a_normed, vec_a0b0)
                vec_a0b1 = b1 - a0
                d1 = np.dot(a_normed, vec_a0b1)
                # Is segment B before A?
                if d0 <= 0 >= d1:
                    if np.absolute(d0) < np.absolute(d1):
                        if len_a0b0 < cutoff:
                            hits.append(potential_interaction)
                            log.debug(
                                "Parallel lines (B first, d0 small) are close: %s", (potential_interaction))
                        else:
                            log.debug(
                                "Parallel lines (B first, d0 small) are far: %s", (potential_interaction))
                        continue
                    if ftuv.magnitude(b1 - a0) < cutoff:
                        hits.append(potential_interaction)
                        log.debug(
                            "Parallel lines (B first, d0 big) are close: %s", (potential_interaction))
                    else:
                        log.debug(
                            "Parallel lines (B first, d0 big) are far: %s", (potential_interaction))
                    continue
                # Is segment B after A?
                elif d0 >= np.asscalar(magnitudes[i]) <= d1:
                    if np.absolute(d0) < np.absolute(d1):
                        if ftuv.magnitude(b0 - a1) < cutoff:
                            hits.append(potential_interaction)
                            log.debug(
                                "Parallel lines (A first, d0 small) are close: %s", (potential_interaction))
                        else:
                            log.debug(
                                "Parallel lines (A first, d0 small) are far: %s", (potential_interaction))
                        continue
                    if ftuv.magnitude(b1 - a1) < cutoff:
                        hits.append(potential_interaction)
                        log.debug(
                            "Parallel lines (A first, d0 big) are close: %s", (potential_interaction))
                    else:
                        log.debug(
                            "Parallel lines (A first, d0 big) are far: %s", (potential_interaction))
                    continue
                if ftuv.magnitude(((d0 * a_normed) + a0) - b0) < cutoff:
                    hits.append(potential_interaction)
                    log.debug("Parallel lines (ELSE) are close: %s",
                              (potential_interaction))
                else:
                    log.debug("Parallel lines (ELSE) are far: %s",
                              potential_interaction)
                continue
            # Lines criss-cross: Calculate the dereminent

            # np.linalg.det([vec_a0b0, b_normed, cross])
            det0 = ftuv.det3x3(np.array([vec_a0b0, b_normed, cross]))
            # np.linalg.det([vec_a0b0, a_normed, cross])
            det1 = ftuv.det3x3(np.array([vec_a0b0, a_normed, cross]))

            t0 = det0 / denom
            t1 = det1 / denom

            pA = a0 + (a_normed * t0)
            pB = b0 + (b_normed * t1)

            # Clamp results to line segments if needed
            pA_ = pA
            pB_ = pB

            if t0 < 0:
                pA_ = a0
            elif t0 > magnitudes[i]:
                pA_ = a1

            if t1 < 0:
                pB_ = b0
            elif t1 > magnitudes[j]:
                pB_ = b1

            dA = ftuv.magnitude(pA - pA_)
            dB = ftuv.magnitude(pB - pB_)
            if dB > dA:
                dot = np.dot(pB_ - a0, a_normed)

                if dot < 0:
                    pA_ = a0
                elif dot > magnitudes[i]:
                    pA_ = a1
                else:
                    pA_ = a0 + a_normed * dot
            elif dB < dA:
                dot = np.dot(pA_ - b0, b_normed)

                if dot < 0:
                    pB_ = b0
                elif dot > magnitudes[j]:
                    pB_ = b1
                else:
                    pB_ = b0 + b_normed * dot

            pA = pA_
            pB = pB_

            d = ftuv.magnitude(pA - pB)
            log.debug("d {}, cutoff {} for {}".format(
                d, cutoff, potential_interaction))
            if d < cutoff:

                hits.append(potential_interaction)
        return hits

    def rmsd_to(self, other):
        # This import is here to avoid circular imports.
        import forgi.threedee.model.similarity as ftms
        if self._elem_names == other._elem_names:
            return ftms.rmsd(self._coordinates,
                             other._coordinates,
                             self.is_centered & other.is_centered)
        else:
            common_keys = set(self._elem_names.keys()) & set(
                other._elem_names.keys())
            if len(common_keys) == len(self._elem_names):
                rev_lookup = list(x for x in sorted(
                    self._elem_names.keys(), key=self._elem_names.__getitem__))
                other_array = np.array(
                    [coord_line for d in rev_lookup for coord_line in [other[d][0], other[d][1]]])
                return ftms.rmsd(self._coordinates,
                                 other_array,
                                 self.is_centered & other.is_centered)
            else:
                common_keys = list(common_keys)
                this_array = np.array(
                    [coord_line for d in common_keys for coord_line in [self[d][0], self[d][1]]])
                other_array = np.array(
                    [coord_line for d in common_keys for coord_line in [other[d][0], other[d][1]]])
                return ftms.rmsd(this_array,
                                 other_array,
                                 False)

    @profile
    def _indices_for(self, elem_name):
        # Avoiding a range is faster if we know that we have only two keys!
        assert self._coords_per_key == 2
        try:
            i = self._elem_names[elem_name]
        except ValueError:
            raise KeyError("Invalid index {}".format(elem_name))
        return [2 * i, 2 * i + 1]
