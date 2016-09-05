#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip) #future package
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)


import sys
import argparse

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.model.descriptors as ftur
import forgi.threedee.model.similarity as ftms

def main(args):

    if len(args.compareTo)==1:
        cg1 = ftmc.CoarseGrainRNA(args.reference[0])
        cg2 = ftmc.CoarseGrainRNA(args.compareTo[0])
        print (ftms.cg_rmsd(cg1, cg2))
    else:
        print ("{:15}\t{:6}\t{:6}\t{:6}".format("filename","RMSD", "dRMSD", "ACC"))
        ref_cg = ftmc.CoarseGrainRNA(args.reference[0])
        reference = ref_cg.get_ordered_virtual_residue_poss()
        acc_calc = ftms.AdjacencyCorrelation(ref_cg)
        for filename in args.compareTo:
            cg=ftmc.CoarseGrainRNA(filename)
            curr_vress=cg.get_ordered_virtual_residue_poss()
            rmsd  = ftur.rmsd(reference, curr_vress)
            drmsd = ftur.drmsd(reference, curr_vress)
            acc   = ftms.mcc(acc_calc.evaluate(cg))
            print ("{:15}\t{:6.3f}\t{:6.3f}\t{:6.3f}".format(filename[-15:], rmsd, drmsd, acc))


parser = argparse.ArgumentParser(description="Calculate the RMSD between two or more coarse grain models")
parser.add_argument('reference', nargs=1, help="The reference *.cg/*.coord file holding the coarse grained model", type=str)
parser.add_argument('compareTo', nargs='+',  help="Other coarse grained models which will be compared to the first one.", type=str)
if __name__ == '__main__':    
    args=parser.parse_args()
    main(args)


