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
    if descriptor not in AVAILABLE_DESCRIPTORS:
        raise ValueError("Descriptor {} not available.".format(descriptor))
    if descriptor == "rog":
        if domain:
            return np.array([ ftmd.radius_of_gyration(cg.get_poss_for_domain(domain, "vres")) for cg in ensemble ])
        else:
            return np.array([ cg.radius_of_gyration() for cg in ensemble ])
    elif descriptor == "anisotropy":
        if domain:
            return np.array([ ftmd.anisotropy(cg.get_poss_for_domain(domain, "vres")) for cg in ensemble ])
        else:
            return np.array([ ftmd.anisotropy(cg.get_ordered_stem_poss()) for cg in ensemble ])
    elif descriptor == "asphericity":
        if domain:
            return np.array([ ftmd.asphericity(cg.get_poss_for_domain(domain, "vres")) for cg in ensemble ])
        else:
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
    return autocorrelate_data(y, mean)

def autocorrelate_data(y, mean=None):
    if mean is None:
        mean = np.mean(y)
    yunbiased = y - mean
    ynorm = np.sum(yunbiased**2)
    corr=np.correlate(yunbiased, yunbiased, mode="same")/ynorm
    return corr[len(y)/2:]

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
        if isinstance(i, int):
            return self._cgs[i]
        if isinstance(i, slice):
            if i.step is not None:
                raise IndexError("Steps are not supperted when using slice notation "
                                 "on Ensemble objects")
            try:
                return EnsembleView(self, i.start, i.stop)
            except ValueError as e:
                raise IndexError(e)
        raise IndexError("Unsupported index {}"/format(i))
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
        self._ensemble   = ens
        self._start = start
        self._end   = end
    def __len__(self):
        return self._end-self._start
    def __getitem__(self, i):
        if i>len(self):
            raise IndexError(i)
        return self._ensemble[self._start+i]
    def get_descriptor(self, descriptor, domain=None):
        if domain is None:
            if descriptor not in self._ensemble._descriptors:
                return _get_descriptor(self, descriptor, domain)
            else:
                return self._ensemble._descriptors[descriptor][self._start:self._end]
        else: #No caching with domains.
            return _get_descriptor(self, descriptor, domain)

def ensemble_from_filenames(filenames):
    cgs=[]
    for fn in filenames:
        cgs.append(ftmc.CoarseGrainRNA(fn))
    return Ensemble(cgs)

if __name__ == "__main__":
    ens = ensemble_from_filenames(sys.argv[1:])
    corrR = autocorrelation(ens, "rog")
    corrA = autocorrelation(ens, "anisotropy")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(4,2)
    x = np.array(list(range(len(corrR)))) 
    print(len(ens.get_descriptor("rog")), " samples")
    rog = ens.get_descriptor("rog")
    anisotropy = ens.get_descriptor("anisotropy")
    #Full Data
    axA = ax[0,0].twinx()
    ax[0,0].plot(list(range(len(ens))), rog, label="ROG", color="blue")    
    axA.plot(list(range(len(ens))), anisotropy, label="Anisotropy", color="red")    
    ax[0,0].set_ylabel("ROG")    
    axA.set_ylabel("anisotropy")
    ax[0,0].set_xlabel("step")
    axA.legend()

    axA2 = ax[0,1].twinx()    
    axA2.plot(x, corrA, color="blue")
    axA2.plot([0,len(corrA)], [0,0], color="black")
    ax[0,1].plot(x, corrR, color="red")
    ax[0,1].set_ylabel("autocorrelation (ROG)")
    ax[0,1].set_xlabel("delta step")
    axA2.set_ylabel("autocorrelation (anisotropy)")

    # Change in full Data
    delta_rog = rog[1:]-rog[:-1]
    delta_anisotropy = rog[1:]-rog[:-1]
    ac_dROG = autocorrelate_data(delta_rog)
    ax[1,0].plot(list(range(len(delta_rog))), delta_rog, label="Detla ROG", color="blue")
    ax[1,1].plot(list(range(len(ac_dROG))), ac_dROG, color="blue")    
    ax[1,1].plot([0,len(ac_dROG)], [0,0], color="black")
    ax[1,0].set_ylabel("Delta ROG")
    ax[1,1].set_ylabel("autocorrelation (Delta ROG)")
    ax[1,0].set_xlabel("step")
    ax[1,1].set_xlabel("lag")
    start = ens[:50]
    rog_start = start.get_descriptor("rog")
    delta_rog_start = rog_start[1:]-rog_start[:-1]
    ac_dROGstart = autocorrelate_data(delta_rog_start)
    ax[1,1].plot(list(range(len(ac_dROGstart))), ac_dROGstart, color="green", label="Delta ROG (first 50 steps)")    
    ax[1,1].plot([0,len(ac_dROG)], [0,0], color="black")
    ax[1,1].legend()
    #ROG
    domains =  ens[0].get_domains()
    for key in domains:
        for i, domain in enumerate(domains[key]):
            if len(domain)>3:
                ax[2,0].plot(np.array(list(range(len(ens)))), ens.get_descriptor("rog", domain), label="{}{}".format(key,i))
                ax[2,1].plot(x, autocorrelation(ens, "rog", domain))
    ax[2,1].plot([0,len(corrR)], [0,0])
    ax[2,0].legend()
    ax[2,0].set_ylabel("ROG for domains")
    ax[2,0].set_xlabel("step")
    ax[2,1].set_ylabel("autocorrelation (ROG)")
    ax[2,1].set_xlabel("delta step")
    #ROG
    for key in domains:
        for i, domain in enumerate(domains[key]):
            if len(domain)>3:
                ax[3,0].plot(np.array(list(range(len(ens)))), ens.get_descriptor("anisotropy", domain), label="{}{}".format(key,i))
                ax[3,1].plot(x, autocorrelation(ens, "anisotropy", domain))
    ax[3,1].plot([0,len(corrR)], [0,0])
    ax[3,0].set_ylabel("anisotropy for domains")
    ax[3,0].set_xlabel("step")
    ax[3,1].set_ylabel("autocorrelation (anisotropy)")
    ax[3,1].set_xlabel("delta step")
    plt.show()
