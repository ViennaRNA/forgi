#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip) #future package
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

from collections import Mapping, defaultdict
import itertools as it
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.manifold import MDS
import time
import forgi.threedee.model.descriptors as ftmd
import scipy.stats
import matplotlib.pyplot as plt

import logging
log = logging.getLogger(__name__)

np_nans = lambda *args, **kwargs: np.ones(*args, **kwargs) * np.nan
    
class Ensemble(Mapping):
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
            
        ############## Caching of some descriptors ############################
        
        #The rmsd matrix. nan means the value needs to be calculated and will then be stored here.
        self._rmsd = np.ones((len(self._cgs), len(self._cgs)))*np.nan        
        for i in range(len(self._rmsd)):
            self._rmsd[i,i]=0.
        # 1D descriptors
        self._descriptors = {}

    def _get_descriptor(self, descr_name):
        if descr_name not in self._descriptors:
            self._descriptors[descr_name] = calculate_descriptor_for(descr_name, 
                                                               self._cgs, 
                                                               *self._get_args_for(descr_name))
            if descr_name == "rmsd_to_last":
                self._rmsd[-1,:]=self._descriptors[descr_name]
                self._rmsd[:, -1]=self._descriptors[descr_name]

        return self._descriptors[descr_name]
    def _add_to_cg_list(self, cg, key):
        """
        During construction of the ensemble, this is used.
        
        Save some bookkeeping variables and store the cg, 
        if it is not identical to the previouse.
        
        :param cg: The coarse grained RNA
        :param key: The index or name under which this cg can be retrieved again.
        """
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
    
    def __iter__(self):
        return iter(self._cg_lookup.keys())
    
    def __len__(self):
        return len(self._cgs)
    
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
        Fill out all empty fields in the rmsd matrix.
        """
        if np.any(np.isnan(self._rmsd)):
            #print(np.where(np.isnan(self._rmsd)))
            log.info("Starting complete rmsd calculation at {}".format(time.time()))
            for i,j in it.combinations(range(len(self)),2):
                if np.isnan(self._rmsd[i,j]):
                    self._rmsd[i,j]=self._rmsd[j,i] = self._cgs[i].coords.rmsd_to(self._cgs[j].coords)
            log.info("Finished complete rmsd calculation at {}".format(time.time()))
        
    def _cluster_dbscan(self):
        """"
        Cluster all structures based on the DBSCAN algorithm 
        using the pairwise RMSD as distance.
        """
        self._calculate_complete_rmsd_matrix()
        db = DBSCAN(eps=np.mean(self._rmsd[0])/3, min_samples=2, metric="precomputed").fit(self._rmsd)
        return db
    
    def _get_args_for(self, descriptor_name):
        """
        Get the arguments that are required to calculate the given descriptor of the ensemble's cgs.
        
        :param descriptor_name: The name of the descriptor to be calculated.
        """
        if descriptor_name == "rmsd_to_reference":
            if self._reference_cg:
                return [self._reference_cg]
            else:
                return [self._cgs[0]]
        elif descriptor_name == "rmsd_to_last":
            return [self._cgs[-1]]
        else:
            return []
        
    def view_delta_rmsd_vs_steps(self):
        self._calculate_complete_rmsd_matrix()
        fig, axes = plt.subplots(2)
        a_rmsd = np_nans(len(self._cg_sequence)//2)
        min_rmsd = np_nans(len(self._cg_sequence)//2)
        max_rmsd = np_nans(len(self._cg_sequence)//2)
        for d in range(len(a_rmsd)):
            l = [ self._rmsd[self._cg_sequence[i], self._cg_sequence[i+d]] for i in range(len(self._cg_sequence)-d) ]
            a_rmsd[d] = sum(l)/len(l)
            min_rmsd[d] = min(l)
            max_rmsd[d] = max(l)    
        for ax in axes:
            ax.set_xlabel("Steps apart")
            ax.set_ylabel("Average RMSD")
            ax.plot(list(range(len(a_rmsd))), a_rmsd, label="Average RMSD")
            ax.plot(list(range(len(min_rmsd))), min_rmsd, label="Minimal RMSD")
            ax.plot(list(range(len(max_rmsd))), max_rmsd, label="Maximal RMSD")
            ax.plot([0, len(max_rmsd)], [np.max(self._rmsd), np.max(self._rmsd)], "-.", label="Maximal RMSD in whole simulation")
            ax.plot([0, len(max_rmsd)], [np.mean(self._rmsd), np.mean(self._rmsd)], "-.", label="Average RMSD in whole simulation")
            ax.legend(prop={'size':6})
        axes[1].set_xlim([0,50])
        
        plt.savefig("rmsd_steps_apart_{}.svg".format(self._cgs[0].name))
        
        plt.clf()
        plt.close()

    def view_2d_embedding(self):
        #http://baoilleach.blogspot.co.at/2014/01/convert-distance-matrix-to-2d.html
        # First cluster all structures based on pairwise RMSD
        db = self._cluster_dbscan()
        labels = db.labels_
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        unique_labels = set(labels)

        #Then calculate the 2D coordinates for our embedding
        mds = MDS(n_components=2, dissimilarity="precomputed", random_state=6)
        results = mds.fit(self._rmsd)
        coords = results.embedding_

        #Now plot
        plt.plot( coords[:,0], coords[:,1], '-', color="blue")

        colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = 'k'
            class_member_mask = (labels == k)

            plt.plot( coords[:,0][class_member_mask & core_samples_mask],
                      coords[:,1][class_member_mask & core_samples_mask],
                    'o', markerfacecolor=col, markeredgecolor='k', markersize=6
                  )
            plt.plot( coords[:,0][class_member_mask & ~core_samples_mask],
                      coords[:,1][class_member_mask & ~core_samples_mask],
                      'o', markerfacecolor=col, markeredgecolor=col, markersize=1
                )
        plt.savefig("embedding_{}.svg".format(self._cgs[0].name))
        plt.clf()
        plt.close()
       


    def view_2d_projection(self, ref_ensemble = None, x = "rmsd_to_reference", y = "rmsd_to_last", cluster = False):
        """
        Plot a 2D projection of the ensemble to the given x and y axis, 
        and visualize the results of clustering with DBSCAN.
        
        :param ref_ensemble: An ensemble or a list of cgs. Plotted as a background in the images.
        :param x: A STRING. The descriptor name used as x axis.
        :param y: A STRING. The descriptor name used as y axis.
        
        Saves the resulting plot as a svg in the current directory.
        """

        # Label the plot
        plt.xlabel(x)
        plt.ylabel(y)

        # First, plot the background (reference ensemble)
        if ref_ensemble is not None:
            log.info("Reference ensemble given")
            ref_x = calculate_descriptor_for(x, ref_ensemble, *self._get_args_for(x))
            ref_y = calculate_descriptor_for(y, ref_ensemble, *self._get_args_for(y))
            
            log.info("Plotting reference")
            plt.plot( ref_x, ref_y, 's', markerfacecolor="green", markeredgecolor='green', markersize=8 )
        else:
            log.info("Reference ensemble missing")
            
        #Get the data for plotting of the ensemble
        data_x = self._get_descriptor(x)
        data_y = self._get_descriptor(y)
        
        if ref_ensemble is not None:
            #Without duplicates
            print("KS-Test without duplicates for {} : {}".format(x, scipy.stats.ks_2samp(data_x, ref_x)))
            print("KS-Test without duplicates for {} : {}".format(y, scipy.stats.ks_2samp(data_y, ref_y)))

            #With duplicates
            full_x = [ data_x[i] for i in self._cg_sequence ] #Correctly account for duplicates
            full_y = [ data_y[i] for i in self._cg_sequence ]
            print("KS-Test for {} : {}".format(x, scipy.stats.ks_2samp(full_x, ref_x)))
            print("KS-Test for {} : {}".format(y, scipy.stats.ks_2samp(full_y, ref_y)))


        if cluster:        
            #In the background, plot lines to show the sampling trajectory
            plt.plot( data_x, data_y, '-', color="blue")

            # Then cluster all structures based on pairwise RMSD
            db = self._cluster_dbscan()
            labels = db.labels_
            core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            unique_labels = set(labels)

            colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
            for k, col in zip(unique_labels, colors):
                if k == -1:
                    # Black used for noise.
                    col = 'k'
                class_member_mask = (labels == k)

                plt.plot( data_x[class_member_mask & core_samples_mask],
                          data_y[class_member_mask & core_samples_mask],
                        'o', markerfacecolor=col, markeredgecolor='k', markersize=6
                      )
                plt.plot( data_x[class_member_mask & ~core_samples_mask],
                      data_y[class_member_mask & ~core_samples_mask],
                      'o', markerfacecolor=col, markeredgecolor=col, markersize=1
                    )
        else:
            #In the background, plot lines to show the sampling trajectory
            plt.plot( data_x, data_y, '-o', color="blue")

        if self._reference_cg:
            plt.plot(calculate_descriptor_for(x, [self._reference_cg], *self._get_args_for(x)),
                     calculate_descriptor_for(y, [self._reference_cg], *self._get_args_for(y)), 
                     "x", color="red", markersize = 12, label="reference" )

        figname = "cluster_{}_{}_{}.svg".format(self._cgs[0].name, x,y)
        plt.savefig(figname)
        log.info("Figure {} created".format(figname))
        plt.clf()
        plt.close()
        
        
class DescriptorCalc(object):
    """
    Helper class to calculate descriptors of Coarse grained RNAs for an ensemble (or a list) of cgs.
    """
    
    @staticmethod
    def rmsd_to_stru(cgs, reference_cg):
        rmsd = np_nans(len(cgs))
        for i, cg in enumerate(cgs):
            rmsd[i] = cg.coords.rmsd_to(reference_cg.coords)
        return rmsd
    @staticmethod
    def rog(cgs):
        rogs = np_nans(len(cgs))
        for i, cg in enumerate(cgs):
            rogs[i] = cg.radius_of_gyration()
        return rogs
    @staticmethod
    def anisotropy(cgs):
        ai = np_nans(len(cgs))
        for i, cg in enumerate(cgs):
            ai[i] = ftmd.anisotropy(cg.get_ordered_stem_poss())
        return ai
    
valid_descriptors = {
    "rmsd_to_reference": DescriptorCalc.rmsd_to_stru,
    "rmsd_to_last": DescriptorCalc.rmsd_to_stru,
    "rog": DescriptorCalc.rog,
    "ROG": DescriptorCalc.rog,
    "anisotropy": DescriptorCalc.anisotropy

    }

    
def calculate_descriptor_for(descriptor_name, cgs, *args):
    """Calculate a descriptor."""
    if descriptor_name not in valid_descriptors:
        raise ValueError("Unknown descriptor {}".format(descriptor_name))
    else:
        return valid_descriptors[descriptor_name](cgs, *args)