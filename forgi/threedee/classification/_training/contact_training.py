from __future__ import print_function, division, absolute_import


"""
Classify unspecific interactions in a very crude way.
Whenever residues are close in the pdb, consider this a contact on the cg-level.

Train a regression model with monotonically decreasing interaction probability
for growing distances.

Don't care about angles.

Train seperate models for different loop types.

First implementation: Use a linear function p=f(dist) which
assumes the maximum of the p vs dist plot at the last maximum of this plot
and 0 at the first bin after the last bin with non-zero probability
(we use p>0.02 instead of p>0 to eliminate outliers).
Considering the noisyness of the data, this seems reasonable
and might be better than least-squares, since we do not want
to take cases, where the probability decreases for shorter
distances close to 0 into account.

"""

import os
import argparse
import warnings
import logging
import itertools as it
import json
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

import Bio.PDB as bpdb

import forgi.utilities.commandline_utils as fuc
import forgi.threedee.utilities.vector as ftuv
import forgi.graph.residue as fgr
log=logging.getLogger(__name__)

################################################################################
### Lets use this as a script
### Run as python -m forgi.threedee.classification._training.contact_training
################################################################################
def main():
    parser = fuc.get_rna_input_parser("Train a classifier for unspecific contacts.",
                                      nargs="+", rna_type="pdb")
    parser.add_argument("--dist", type=float, required=True,
                        help="The cutoff-distance between atoms to classify something as an interaction")
    parser.add_argument("--plot", action="store_true",
                        help="Show plots")
    args = parser.parse_args()
    cgs = fuc.cgs_from_args(args, "+", rna_type="pdb", skip_errors=True)
    all_dists=defaultdict(list)
    all_labels=defaultdict(list)
    for cg in cgs:
        interactions = get_all_interactions(cg, args.dist)
        for interaction in ["ii", "hi", "hh"]:
            dists, label = pdb_contacts(cg, interaction, interactions)
            all_dists[interaction].extend(dists)
            all_labels[interaction].extend(label)
    with open("contact_training_result.json", 'w') as f:
        json.dump({"dists":all_dists, "labels":all_labels}, f)
    for interaction in ["ii", "hi", "hh"]:
        dists=np.array(all_dists[interaction])
        labels=np.array(all_labels[interaction])
        probs = []
        step=5
        for b in np.arange(5,60,step):
            mask=np.where((b-step<=dists)&(dists<b))
            if len(labels[mask]):
                p_interaction = sum(labels[mask])/len(labels[mask])
                print("{}-{}: {}".format(b-step, b, p_interaction))
            else:
                p_interaction=0
            probs.append(p_interaction)
        first_i = 0
        last_i=0
        for i,p in enumerate(probs):
            if p==max(probs):
                first_i=i
            if p>0.02:
                last_i=i
        last_i+=1
        k=-max(probs)/((last_i-first_i-0.5)*step) # From middle of bin to start of bin
        #k*(first_i*step+0.5*step)+d==max(probs)
        d=max(probs)-k*(first_i*step+step*0.5)
        f=lambda x: k*x+d
        x=[]
        for a in np.arange(0, 55, step):
            x.append(a)
            x.append(a+step)
        probs2=[]
        for p in probs:
            probs2.append(p)
            probs2.append(p)
        if args.plot:
            import matplotlib.pyplot as plt
            plt.title(interaction)
            plt.plot(x, probs2)
            plt.plot(np.linspace(0,60,10), f(np.linspace(0,60,10)))
            plt.ylim([0,1])
            plt.show()

def pdb_contacts(cg, interaction_type, all_interactions):
    data  = []
    labels=[]
    for elem1, elem2 in it.combinations(cg.defines, r=2):
        if elem1[0] != interaction_type[0]:
            continue
        if elem2[0]!=interaction_type[1]:
            continue
        if elem1==elem2:
            continue
        if elem1 in cg.incomplete_elements or elem1 in cg.interacting_residues:
            continue
        if elem2 in cg.incomplete_elements or elem1 in cg.interacting_residues:
            continue
        cg_dist=ftuv.vec_distance(*ftuv.line_segment_distance(cg.coords[elem1][0], cg.coords[elem1][1], *cg.coords[elem2]))
        i_pair = tuple(sorted([elem1, elem2]))
        is_interaction = i_pair in all_interactions
        data.append(cg_dist)
        labels.append(int(is_interaction))
    return np.array(data), np.array(labels)

def get_all_interactions(cg, pdb_dist_cutoff):
    atoms = [ a for c in cg.chains.values() for a in c.get_atoms() if a.name[0] in ["C", "N", "S"]]
    kdtree = bpdb.NeighborSearch(atoms)
    pairs = kdtree.search_all(pdb_dist_cutoff, "A")
    res_pair_list=set()
    for a1, a2 in pairs:
        p1 = a1.get_parent()
        p2 = a2.get_parent()
        if p1.id == p2.id:
            continue
        elif p1 < p2:
            res_pair_list.add((p1, p2))
        else:
            res_pair_list.add((p2, p1))
    interacting_elements = set()
    for res1, res2 in res_pair_list:
        try:
            elem1 = cg.get_elem(fgr.RESID(res1.parent.id, res1.id))
            elem2 = cg.get_elem(fgr.RESID(res2.parent.id, res2.id))
        except (ValueError, LookupError):
            log.exception(cg.name)
            continue
        if elem1 != elem2:
            interacting_elements.add(tuple(sorted([elem1, elem2])))
    return interacting_elements

if __name__=="__main__":
    main()
