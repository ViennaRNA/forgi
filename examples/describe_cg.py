# coding: utf-8
from __future__ import print_function, absolute_import, division, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip) 



import argparse
from collections import defaultdict
import forgi.threedee.model.descriptors as ftmd
import forgi.threedee.model.coarse_grain as ftmc
import pandas as pd
import logging
log = logging.getLogger(__name__)
def generateParser():
    parser=argparse.ArgumentParser( description="Collect data about a list of pdb files and store it as a csv.")
    parser.add_argument("files", type=str, nargs="+", help="One or more cg files.")
    parser.add_argument("--verbose", "-v", action="store_true", help="Be verbose.")
    parser.add_argument("--csv", type=str, help="Store dataframe under this filename. (Prints to stdout if not given)")

    return parser


parser = generateParser()
if __name__=="__main__":
    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    data = defaultdict(list)
    for i, filename in enumerate(args.files):
        cg = ftmc.CoarseGrainRNA(filename)
        log.debug("File {} read".format(cg.name))
        data["pdb_name"].append(cg.name)
        data["path"].append(filename)
        data["nt_length"].append(cg.seq_length)
        data["num_cg_elems"].append(len(cg.defines))
        for letter in "smifth":
            data["num_"+letter].append(len([x for x in cg.defines if x[0] ==letter]))
        multiloops = cg.find_mlonly_multiloops()
        descriptors = []
        for ml in multiloops:
            descriptors.append(cg.describe_multiloop(ml))
        #print (descriptors)
        data["open_mls"].append(len([d for d in descriptors if "open" in d]))
        #print(data["open_mls"][-1])
        data["pseudoknots"].append(len([d for d in descriptors if "pseudoknot" in d]))
        data["regular_mls"].append(len([d for d in descriptors if "regular_multiloop" in d]))
        data["total_mls"].append(len(multiloops))
        try:
            data["longest_ml"].append(max(len(x) for x in multiloops))
        except ValueError:
            data["longest_ml"].append(0)
        try:
            data["rog_fast"].append(cg.radius_of_gyration("fast"))
        except:
            data["rog_fast"].append(float("nan"))
            data["rog_vres"].append(float("nan"))
            data["anisotropy_fast"].append(float("nan"))
            data["anisotropy_vres"].append(float("nan"))
            data["asphericity_fast"].append(float("nan"))
            data["asphericity_vres"].append(float("nan"))
        else:
            data["rog_vres"].append(cg.radius_of_gyration("vres"))
            data["anisotropy_fast"].append(ftmd.anisotropy(cg.get_ordered_stem_poss()))
            data["anisotropy_vres"].append(ftmd.anisotropy(cg.get_ordered_virtual_residue_poss()))
            data["asphericity_fast"].append(ftmd.asphericity(cg.get_ordered_stem_poss()))
            data["asphericity_vres"].append(ftmd.asphericity(cg.get_ordered_virtual_residue_poss()))
    df = pd.DataFrame(data)
    if args.csv:
        df.to_csv(args.csv, mode='a')
    else:
        print(df)
