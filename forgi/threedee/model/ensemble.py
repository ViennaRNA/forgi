#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip) #future package
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

from collections import Mapping
import itertools as it
import numpy as np
from sklearn.cluster import DBSCAN
import time
import logging
log = logging.getLogger(__name__)

class Ensemble(object):
    #################### INITIALIZATION
    def __init__(self, cgs, reference_cg = None, sort_key = None):
        """
        An Ensemble is a sequence of Coarse grained RNAs, all of which must correspond 
                    to the same RNA 2D structure.

        :param cgs: An iterable of coarse grain RNAs or a mapping key->cg, 
                    from which the ensemble is constructed.
                    The ensemble may keep references to some of the cgs and 
                    modify them by centering them on the centroid.
                    The ensemble is not guaranteed to keep/copy/reference any of the 
                    cg's properties except the coords and twists.
        :param sort_key: Optional. A function that takes a cg (if cgs is a sequence)
                         or a key (if cgs is a mapping)
        """
        self._reference_cg = reference_cg

        self._cgs=[]
        self._cg_lookup={}
        #The order of cgs, as specified by sort_key
        self._cg_sequence = []
        self._cg_rev_lookup={}
        if isinstance(cgs, Mapping):
            for key in sorted(cgs.keys(), key=sort_key):
                self._add_to_cg_list(cgs[key], key)
        else:
            if sort_key is not None:
                ids = list(sorted(range(len(cgs)), key=lambda x: sort_key(cgs[x])))
            else:
                ids = list(range(len(cgs)))

            for key in ids:            
                self._add_to_cg_list(cgs[key], key)
        #Center all cg coords on their centroid
        for cg in self._cgs:            
            cg.coords.center()
        if self._reference_cg is not None:
            self._reference_cg.coords.center()
        #The rmsd matrix. nan means the value needs to be calculated and will then be stored here.
        self._rmsd = np.ones((len(self._cgs), len(self._cgs)))*np.nan        
        for i in range(len(self._rmsd)):
            self._rmsd[i,i]=0.
        # The rmsd to the reference
        self._ref_rmsd = np.ones(len(self._cgs))*np.nan
    def _add_to_cg_list(self, cg, key):
        log.debug("Adding cg {} with key {}".format(cg.name, key))
        #In the future, we might have to check the rmsd in addition, 
        # if the MCMC will allow rotation/translation
        if not self._cgs or cg.coords != self._cgs[-1].coords or cg.twists != self._cgs[-1].twists:
            if self._cgs:
                log.debug("Coords or twists different!")
            else:
                log.debug("First cg added")
            self._cgs.append(cg)
        else:
            log.debug("It is the same as the previous")
        #Else: only store lookup entry pointing to previous cg
        self._cg_lookup[key]=len(self._cgs)-1
        self._cg_rev_lookup[len(self._cgs)-1]=key
        self._cg_sequence.append(len(self._cgs)-1)
    ######################## MEMBER LOOKUP    
    def __getitem__(self, key):
        """Using a key like used upon ensemble construction, retrieve the corresponding element"""
        return self._cgs[self._cg_lookup[key]]
    def at_timestep(self, timestep):
        """
        Access an ensemble member by the timestep.
             
        In contrast to __getitem__, which uses the index supplied upon ensemble creation,
        this uses the sort order specified by sort_key in the __init__ function.
        It is useful if the ensemble is seen as a trajectory. 
        This method gives the ith frame of this trajectory.
        
        :param timestep: The number of the frame (cg) in the trajectory that should be retrieved.
        :returns: A coarse-grained RNA.
        """
        return self._cgs[self._cg_sequence[timestep]]
    ####################### RMSD and RMSD based calculations
    def rmsd_between(self, key1, key2, mode="key"):
        """
        Return (and cache) the rmsd between two structures.
        
        :param key1, key2: Two keys to reference two cgs
        :param mode: "key" or "timestep". Whether the keys are timesteps or the keys 
                     used by __getitem__
        :returns: the rmsd as float
        """
        if mode=="key":
            i = self._cg_lookup[key1]
            j = self._cg_lookup[key2]
        elif mode == "timestep":
            i = self._cg_sequence[key1] 
            j = self._cg_sequence[key2]
        else: 
            raise ValueError("Invalid mode {}".format(mode))
        if np.isnan(self._rmsd[i,j]):
                self._rmsd[i,j]=self._rmsd[j,i] = self._cgs[i].coords.rmsd_to(self._cgs[j].coords)
        return self._rmsd[i,j]
    def _calculate_complete_rmsd_matrix(self):
        """
        Fill out all empty fields in the rmsd matrix
        """
        for i,j in it.combinations(range(len(self._cgs,2))):
            if np.isnan(self._rmsd[i,j]):
                self._rmsd[i,j]=self._rmsd[j,i] = self._cgs[i].coords.rmsd_to(self._cgs[j].coords)
    def _calculate_complete_ref_rmsd(self):
        log.info("Starting complete rmsd calculation at {}".format(time.time()))
        for i in range(len(self._cgs)):
            self._ref_rmsd[i]=self._cgs[i].coords.rmsd_to(self._reference_cg.coords)
        log.info("Finished complete rmsd calculation at {}".format(time.time()))
    def _cluster_dbscan(self):
        self._calculate_complete_rmsd_matrix()
        db = DBSCAN(eps=1, min_samples=10, metric="precomputed", n_jobs = 4).fit(self._rmsd)
        return db
        
    def view_db_clustering(self):
        import matplotlib.pyplot as plt
        db = self._cluser_dbscan()
        labels = db.labels_
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
    
        unique_labels = set(labels)
        colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
        if self._reference_cg is not None:
            self._calculate_complete_ref_rmsd()
            ref = self._ref_rmsd        
            plt.xlabel("RMSD to reference")
        else:
            ref = self._rmsd[-1]
            plt.xlabel("RMSD to last structure")
            
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = 'k'
            class_member_mask = (labels == k)
            
            plt.plot( ref[class_member_mask & core_samples_mask],
                      self._rmsd[0][class_member_mask & core_samples_mask],
                    'o', markerfacecolor=col, markeredgecolor='k', markersize=8
                  )
            plt.plot( ref[class_member_mask & ~core_samples_mask],
                  self._rmsd[0][class_member_mask & ~core_samples_mask],
                  'o', markerfacecolor=col, markeredgecolor='k', markersize=4
                )

        plt.ylabel("RMSD to start")
        
        
        
        
        
        
