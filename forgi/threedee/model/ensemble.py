from __future__ import absolute_import, unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip)


from collections import MutableSequence, Sequence
import sys
import numpy as np
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.descriptors as ftmd

AVAILABLE_DESCRIPTORS = ["rog", "anisotropy", "asphericity"]

def _get_descriptor(ensemble, descriptor, domain=None):
    """
    :param ensemble: an Ensemble or EnsembleView object
    :param descriptor: A STRING. One of AVAILABLE_DESCRIPTORS
    :param domain: An iterable of cg element names or None (whole cg)
    :returns: A np.array
    """
    if measure not in AVAILABLE_DESCRIPTORS:
        raise ValueError("Descriptor {} not available.".format(descriptor))
    if domain != None: 
        raise NotImplementedError("Domain not yet implemented for get_descriptor")
    if descriptor == "rog":
        return np.array([ cg.radius_of_gyration() for cg in ensemble ])
    elif descriptor == "anisotropy":
        return np.array([ ftmd.anisotropy(cg.get_ordered_stem_poss()) for cg in ensemble ])
    elif descriptor == "asphericity":
        return np.array([ ftmd.asphericity(cg.get_ordered_stem_poss()) for cg in ensemble ])

def autocorrelation(ensemble, descriptor="rog", domain = None, mean = None):
    """
    Return the normalized autocorrelation as a 1D array for the given measure along the
    sequence of structures in this ensemble.
    
    :param ensemble: A ensemble or EnsembleView object
    :param descriptor: A STRING. One of AVAILABLE_DESCRIPTORS
    :param domain: An iterable of cg element names or None (whole cg)
    :param mean: A FLOAT or None. 
                  If this is None, all datapoints will be shifted by their mean 
                  before calculating the autocorrelation.
                  If this is numeric, datapoints will be shifted by this value instead.
    #http://stackoverflow.com/a/17090200/5069869
    """
    y = ensemble.get_descriptor(descriptor, domain)
    if mean is None:
        mean = np.mean(y)
    yunbiased = y - mean
    ynorm = np.sum(yunbiased**2)
    corr=np.correlate(yunbiased, yunbiased, mode="same")/ynorm
    return corr[len(self.rogs)/2:]

class Ensemble(Sequence):
    def __init__(self, cgs):
        """
        :param cgs: An iterable of coarse grain RNAs, all of which must correspond 
                    to the same RNA 2D structure. They should be ordered.
        """
        self._cgs = cgs
        #Cached Data
        self._descriptors = {}

    # Methods for accessing the stored cg structures
    def __getitem__(self, i):
        if len(i)>1:
            raise IndexError("Only a single index/ slice is allowed")
        if isinstance(i, int):
            return self._cgs[i]
        elif isinstance(i, slice_):
            if slice_.step is not None:
                raise IndexError("Steps are not supperted when using slice notation on Ensemble objects")
            try:
                return EnsembleView(self, slice_.start, slice_.stop)
            except ValueError as e:
                raise IndexError(e)
    def __len__(self):
        return len(self._cgs)
    """# Methods for manipulating the list of cgs
    def __setitem__(self, i, cg):
        self._cgs[i]=cg
        self.update()
    def __delitem__(self, i):
        del self._cgs[i]
        self.update()
    def insert(self, i, cg):
        self._cgs.insert(i,cg)
        self.update()
    #"""
    def update(self):
        """
        Clear all cached values.
        
        Has to be called, whenever a CoarseGrainRNA that is 
        part of the ensemble is modified in palce.

        ..note::
            This is called automatically, when CoarseGrainRNAs are added to/ 
            removed from the ensemble.
        """
        self._descriptors = {}

    def get_descriptor(self, descriptor, domain=None):
        if domain is None:
            if descriptor not in self._descriptors:
                self._descriptors[descriptor] = _get_descriptor(self, descriptor, domain)
            return self._descriptors[descriptor]
        else: #No caching with domains.
            return _get_descriptor(self, descriptor, domain)

class EnsembleView(Sequence):
    def __init__(self, ens, start, end):
        if start is None:
            start = 0
        if end is None:
            end = len(ens)
        if not isinstance(start, int):
            raise ValueError("Start of EnsembleView must be integer.")
        if not isinstance(end, int):
            raise ValueError("End of EnsembleView must be integer.")
        if start<0 or start>end or end<0:
            raise ValueError("Start and End of slice EnsembleView must be positive and end>start!.")
        if end>len(ens):
            end=len(ens)
        self.ensemble   = ens
        self.start = start
        self.end   = end
    def __len__(self):
        return self.end-self.start
    def __getitem__(self, i):
        if i>len(self):
            raise IndexError(i)
        return self.ensemble[start+i]
    def get_descriptor(self, descriptor, domain=None):
        if domain is None:
            if descriptor not in self.ens._descriptors:
                return _get_descriptor(self, descriptor, domain)
            else:
                return self.ens._descriptors[descriptor]
        else: #No caching with domains.
            return _get_descriptor(self, descriptor, domain)

def ensemble_from_filenames(filenames):
    cgs=[]
    for fn in filenames:
        cgs.append(ftmc.CoarseGrainRNA(fn))
    return Ensemble(cgs)

if __name__ == "__main__":
    ens = ensemble_from_filenames(sys.argv[1:])
    corrR = ens.autocorrelation()
    corrA = ens.autocorrelation("anisotropy")

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2,2)
    x = np.array(list(range(len(corrR))))

    #ROG
    ax[0,0].plot(np.array(list(range(len(ens.rogs)))), ens.rogs)
    ax[0,1].plot(x, corrR)
    ax[0,0].set_ylabel("ROG")
    ax[0,0].set_xlabel("step")
    ax[0,1].set_ylabel("autocorrelation (ROG)")
    ax[0,1].set_xlabel("delta step")
    #Anisotropy
    ax[1,0].plot(np.array(list(range(len(ens.anisotropies)))), ens.anisotropies)
    ax[1,1].plot(x, corrA)
    ax[1,0].set_ylabel("anisotropy")
    ax[1,0].set_xlabel("step")
    ax[1,1].set_ylabel("autocorrelation (anisotropy)")
    ax[1,1].set_xlabel("delta step")
    plt.show()
