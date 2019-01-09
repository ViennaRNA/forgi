#!/usr/bin/env python
from __future__ import print_function

import logging
import os.path
import sys
import numpy as np

import forgi.threedee.model.similarity as ftms
import forgi.utilities.commandline_utils as fuc
import forgi.threedee.utilities.pdb as ftup

log = logging.getLogger(__name__)


def get_parser():
    parser = fuc.get_rna_input_parser("Compare two RNA 3d structures based on RMSD, ACC"
                                      " or other measures", 2, "3d", enable_logging=True)
    output_group = parser.add_argument_group("Controlling output",
                                             description="Control, based on which measure the two "
                                             "structures will be compared")
    output_group.add_argument("--acc", help="Compare based on the Adjacency correlation coefficient"
                                            " ACC", action="store_true")
    output_group.add_argument(
        "--rmsd", help="Compare based on CG-RMSD", action="store_true")
    output_group.add_argument(
        "--pdb-rmsd", help="Compare based on PDB-RMSD", action="store_true")
    return parser


def main(args):
    with fuc.hide_traceback():
        cg1, cg2 = fuc.cgs_from_args(
            args, rna_type="3d", enable_logging=True)

        if not (args.acc or args.rmsd or args.pdb_rmsd):
            showall = True
        else:
            showall = False
        if showall or args.acc:
            if cg1.defines != cg2.defines:
                if args.acc:
                    print(
                       "Cannot compare two 3d structures that do not correspond to the same RNA.")
                    sys.exit(1)
            else:
                adj = ftms.AdjacencyCorrelation(cg1)
                print("ACC:\t{:.3f}".format(ftms.mcc(adj.evaluate(cg2))))
        if showall or args.rmsd:
            print("RMSD:\t{:.3f}".format(ftms.cg_rmsd(cg1, cg2)))
        if showall or args.pdb_rmsd:
            if not pdb_rmsd(cg1, cg2):
                # If --pdb-rmsd was not given, just don't print it.
                # If it was given, we exit with non-zero exit status.
                if args.pdb_rmsd:
                    print(
                        "Cannot calculate PDB-RMSD: The two files do not contain the same chains.")
                    sys.exit(1)

def pdb_rmsd(cg1, cg2):
    common_chains = set(cg1.chains.keys()) & set(cg2.chains.keys())
    reslist1 = []
    reslist2 = []
    for chain in common_chains:
        reslist1.extend(cr for cr in cg1.chains[chain].get_list()
                         if cr.resname.strip() in ftup.RNA_RESIDUES)
        reslist2.extend(cr for cr in cg2.chains[chain].get_list()
                         if cr.resname.strip() in ftup.RNA_RESIDUES)
    if common_chains:
        print("PDB-RMSD (chains {}):\t{:.3f}".format("-".join(common_chains),
              ftup.residuelist_rmsd(reslist1, reslist2, sidechains=True)[1]))
        return True
    return False


parser = get_parser()
if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
