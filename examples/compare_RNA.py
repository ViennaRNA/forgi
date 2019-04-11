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
    parser = fuc.get_rna_input_parser("Compare RNA 3d structures based on RMSD, ACC"
                                      " or other measures. If more than 2 filenames are given, compare all files to the first.", "+", "3d", enable_logging=True)
    output_group = parser.add_argument_group("Controlling output",
                                             description="Control, based on which measure the two "
                                             "structures will be compared")
    output_group.add_argument("--acc", help="Compare based on the Adjacency correlation coefficient"
                                            " ACC", action="store_true")
    output_group.add_argument(
        "--rmsd", help="Compare based on CG-RMSD", action="store_true")
    output_group.add_argument(
        "--stem-rmsd", help="Compare based on coarse grained stem RMSD", action="store_true")
    output_group.add_argument(
        "--pdb-rmsd", help="Compare based on PDB-RMSD", action="store_true")
    return parser


def main(args):
    with fuc.hide_traceback():
        cgs = fuc.cgs_from_args(
            args, rna_type="3d", enable_logging=True)

        if len(cgs)<2:
            raise ValueError("At least 2 RNA structures required for comparison")
        cg1=cgs[0]
        if not (args.acc or args.rmsd or args.pdb_rmsd or args.stem_rmsd):
            showall = True
        else:
            showall = False
        for cg2 in cgs[1:]:
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
            if showall or args.stem_rmsd:
                try:
                    print("Stem:\t{:.3f}".format(ftms.cg_stem_rmsd(cg1, cg2)))
                except ftms.Incompareable:
                    if args.stem_rmsd:
                        raise
                    # Otherwise do nothing
            if showall or args.pdb_rmsd:
                if not pdb_rmsd(cg1, cg2):
                    # If --pdb-rmsd was not given, just don't print it.
                    # If it was given, we exit with non-zero exit status.
                    if args.pdb_rmsd:
                        print("Cannot calculate PDB-RMSD: "
                              "The two files do not contain the same chains.")

def pdb_rmsd(cg1, cg2):
    if len(cg1.chains)==1==len(cg2.chains):
        reslist1 = []
        reslist2 = []
        chain1, = cg1.chains.values()
        chain2, = cg2.chains.values()
        resnames1=[ r.resname.strip() for r in chain1 ]
        resnames2=[ r.resname.strip() for r in chain2 ]
        log.info("Resname lists are %s, %s", resnames1, resnames2)
        if resnames1 == resnames2:
            for r in chain1:
                if r.resname.strip() in ftup.RNA_RESIDUES:
                    reslist1.append(r)
            for r in chain2:
                if r.resname.strip() in ftup.RNA_RESIDUES:
                    reslist2.append(r)
            print("PDB-RMSD (single chain):\t{:.3f}".format(
                      ftup.residuelist_rmsd(reslist1, reslist2, sidechains=True)[1]))
            return True

        else:
            log.info("Resname lists are different")
            reslist1, reslist2 = common_reslists(chain1, chain2)
            print("PDB-RMSD (single chain):\t{:.3f}".format(
                      ftup.residuelist_rmsd(reslist1, reslist2, sidechains=True)[1]))
            return True
    else:
        log.info("Chain lengths are %s and %s",len(cg1.chains), len(cg2.chains) )
        common_chains = set(cg1.chains.keys()) & set(cg2.chains.keys())
        reslist1 = []
        reslist2 = []
        for chain in common_chains:
            ext1, ext2 = common_reslists(cg1.chains[chain], cg2.chains[chain])
            reslist1.extend(ext1)
            reslist2.extend(ext2)
        if common_chains:
            print("PDB-RMSD (chains {}):\t{:.3f}".format("-".join(common_chains),
                  ftup.residuelist_rmsd(reslist1, reslist2, sidechains=True)[1]))
            return True
        return False


def common_reslists(chain1, chain2):
    reslist1 = []
    reslist2 = []
    start1 = chain1.get_list()[0].id[1]
    start2 = chain2.get_list()[0].id[1]
    offset=start2-start1
    mismatches=float("inf")
    for offs in [0, offset]:
        rl1, rl2, m = _common_reslist_with_offset(chain1, chain2, offs)
        if m<mismatches:
            mismatches=m
            reslist1, reslist2 = rl1, rl2
    if mismatches:
        log.warning("There were %s point mutations between the two chains", mismatches)
    if mismatches>8:
        raise ValueError("Could not find mapping between chains!")
    return reslist1, reslist2

def _common_reslist_with_offset(chain1, chain2, offset):
    log.info("Trying mapping with offset %s", offset)
    reslist1 = []
    reslist2 = []
    for cr2 in chain2.get_list():
        mapped_id = cr2.id[0], cr2.id[1]-offset, cr2.id[2]
        if mapped_id not in chain1:
            log.info("Residue %s (%s) [->%s] not in chain1", cr2.id, cr2.resname.strip(), mapped_id)
            log.info("Res %s has atoms %s", mapped_id, cr2.get_list())
    mismatch=0
    for cr1 in chain1.get_list():
        mapped_id = cr1.id[0], cr1.id[1]+offset, cr1.id[2]
        cr1_name = cr1.resname.strip()
        try:
            cr2 = chain2[mapped_id]
        except KeyError as e:
            log.info("Residue %s (%s) [->%s] not in chain2", cr1.id, cr1_name, mapped_id)
            continue
        cr2_name = cr2.resname.strip()
        if cr1_name != cr2_name:
            mismatch+=1
            log.info("Residue %s [->%s](chains %s, %s) has different name in the "
                        "two chains: %s!=%s", cr1.id, mapped_id, chain1.id, chain2.id, cr1_name, cr2_name)
        if cr1_name in ftup.RNA_RESIDUES and cr2_name in ftup.RNA_RESIDUES:
            reslist1.append(cr1)
            reslist2.append(cr2)

    return reslist1, reslist2, mismatch



parser = get_parser()
if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
