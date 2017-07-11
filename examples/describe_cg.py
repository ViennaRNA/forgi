# coding: utf-8
from __future__ import print_function, absolute_import, division, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip)



import argparse
from collections import defaultdict
import forgi.utilities.commandline_utils as fuc
import forgi.threedee.model.descriptors as ftmd
import forgi.threedee.model.coarse_grain as ftmc
import pandas as pd
import logging
log = logging.getLogger(__name__)

def generateParser():
    parser=fuc.get_rna_input_parser("Collect data about a list of rna files and store it as a csv.",
                                     nargs="+", rna_type="any", enable_logging=True)
    parser.add_argument("--csv", type=str, help="Store dataframe under this filename. (Prints to stdout if not given)")
    parser.add_argument("-k", "--keys", type=str, help="Only print the following properties. "
                                        "(A comma-seperated list of column headers, e.g. rog_vres")

    return parser


parser = generateParser()
if __name__=="__main__":
    args = parser.parse_args()
    cgs = fuc.cgs_from_args(args, "+", "any", enable_logging=True)
    data = defaultdict(list)
    for i, cg in enumerate(cgs):
        data["name"].append(cg.name)
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
    if args.keys:
        allowed_keys = args.keys.split(",")
        for key in list(data.keys()):
            if key not in allowed_keys:
                del data[key]
    df = pd.DataFrame(data)
    if args.csv:
        df.to_csv(args.csv, mode='a')
    else:
        print(df)
