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

import numpy as np
import matplotlib.pyplot as plt


import forgi.utilities.commandline_utils as fuc
import forgi.threedee.utilities.vector as ftuv

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

    args = parser.parse_args()
    cgs = fuc.cgs_from_args(args, "+", rna_type="pdb")
    for interaction in ["ii", "hi", "hh"]:
        all_dists, all_labels=[],[]
        for cg in cgs:
            dists, label = pdb_contacts(cg, interaction, args.dist)
            all_dists.extend(dists)
            all_labels.extend(label)
            step=5
        all_dists=np.array(all_dists)
        all_labels=np.array(all_labels)
        probs = []
        for b in np.arange(5,60,step):
            mask=np.where((b-step<=all_dists)&(all_dists<b))
            if len(all_labels[mask]):
                p_interaction = sum(all_labels[mask])/len(all_labels[mask])
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

        plt.plot(x, probs2)
        plt.plot(np.linspace(0,60,10), f(np.linspace(0,60,10)))
        plt.ylim([0,1])
        plt.show()
def pdb_contacts(cg, interaction_type, pdb_dist_cutoff=3):
    data  = []
    labels=[]
    for elem1, elem2 in it.combinations(cg.defines, r=2):
        if elem1[0] != interaction_type[0]:
            continue
        if elem2[0]!=interaction_type[1]:
            continue
        if elem1==elem2:
            continue
        cg_dist=ftuv.vec_distance(*ftuv.line_segment_distance(cg.coords[elem1][0], cg.coords[elem1][1], *cg.coords[elem2]))
        is_interaction = are_elements_interacting(cg, elem1, elem2, pdb_dist_cutoff)
        data.append(cg_dist)
        labels.append(is_interaction)
    return np.array(data), np.array(labels)

def are_elements_interacting(cg, elem1, elem2, pdb_dist_cutoff):
    # Inefficient, but who cares
    for r1, r2 in it.product(cg.define_residue_num_iterator(elem1, seq_ids=True),
                             cg.define_residue_num_iterator(elem2, seq_ids=True)):
        res1 = cg.chains[r1.chain][r1.resid]
        res2 = cg.chains[r2.chain][r2.resid]
        for atom1 in res1:
            for atom2 in res2:
                if ftuv.vec_distance(atom1.get_coord(), atom2.get_coord())<pdb_dist_cutoff:
                    return True
    return False

if __name__=="__main__":
    main()
