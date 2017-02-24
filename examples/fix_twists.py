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
import forgi.threedee.utilities.vector as ftuv


parser = argparse.ArgumentParser(description="Read a CoarseGrainRNA-File that may contain twists not normal to the stem axis, and put out a corrected version of this cg file.")
parser.add_argument('rna', help="The *.cg/*.coord file holding the coarse grained model", type=str)


def main(args):
    cg = ftmc.CoarseGrainRNA(args.rna)
    for stem in cg.stem_iterator():            
        stem_vec = cg.coords.get_direction(stem)
        twist_vec = cg.twists[stem][0]
        try:
            ftuv.create_orthonormal_basis(stem_vec, twist_vec)
        except ValueError:
            cg.twists[stem] = ftuv.vector_rejection(twist_vec, stem_vec), cg.twists[stem][1]
        twist_vec = cg.twists[stem][1]
        try:
            ftuv.create_orthonormal_basis(stem_vec, twist_vec)
        except ValueError:
            cg.twists[stem] = cg.twists[stem][0], ftuv.vector_rejection(twist_vec, stem_vec)
    try:
        cg.add_all_virtual_residues()
    except:
        assert False
    print(cg.to_cg_string())
    
if __name__ == '__main__':    
    args=parser.parse_args()
    main(args)


