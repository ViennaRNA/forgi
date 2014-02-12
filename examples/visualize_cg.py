#!/usr/bin/python

import numpy as np
import sys
import subprocess as sp
import tempfile as tf

import forgi.threedee.model.coarse_grain as cmg
import forgi.utilities.debug as cud
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.pdb as cup
import forgi.threedee.utilities.rmsd as ftur
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.visual.pymol as cvp

from optparse import OptionParser

def align_cgs(cgs):
    '''
    Align each coarse grain RNA to the first one.
    
    @param cgs: A list of CoarseGrainRNA structures.
    @return: Nothing, the cgs are modified in place
    '''
    centroid0 = ftuv.get_vector_centroid(ftug.bg_virtual_residues(cgs[0]))
    crds0 = ftuv.center_on_centroid(ftug.bg_virtual_residues(cgs[0]))

    for cg in cgs:
        centroid1 = ftuv.get_vector_centroid(ftug.bg_virtual_residues(cg))
        crds1 = ftuv.center_on_centroid(ftug.bg_virtual_residues(cg))

        rot_mat = ftur.optimal_superposition(crds0, crds1)
        for k in cg.coords.keys():
            cg.coords[k] = (np.dot(rot_mat, cg.coords[k][0] - centroid1),
                            np.dot(rot_mat, cg.coords[k][1] - centroid1))

            if k[0] == 's':
                cg.twists[k] = (np.dot(rot_mat, cg.twists[k][0]),
                                np.dot(rot_mat, cg.twists[k][1]))


def main():
    usage = """
    ./visualize_cg.py cg_file

    Display the coarse-grain representation of a structure in pymol.
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')
    parser.add_option('-r', '--longrange', dest='longrange', default=False, action='store_true', help="Display long-range interactions")
    parser.add_option('-l', '--loops', dest='loops', default=True, action='store_false', help="Don't display the coarse-grain hairpin loops")
    parser.add_option('-c', '--cones', dest='cones', default=False, action='store_true', help="Display cones that portrude from the stems")
    parser.add_option('-x', '--text', dest='text', default=False, action='store_true', help="Add labels to the figure.")
    parser.add_option('-a', '--align', dest='align', default=False, action='store_true', help='Align all of the structures with the first')
    parser.add_option('-e', '--encompassing-stems', dest='encompassing_stems', default=False, action='store_true', help='Show the big stems that encompass the colinear ones.')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    pp = cvp.PymolPrinter()
    pp.add_loops = options.loops
    pp.draw_cones = options.cones
    #sys.exit(1)
    pp.add_longrange = options.longrange
    pp.print_text = options.text
    pp.encompassing_stems = options.encompassing_stems

    cgs = []
    for a in args:
        cgs += [cmg.from_file(a)]

    if options.align:
        align_cgs(cgs)

    for i,cg in enumerate(cgs):
        if i > 0:
            pp.override_color = 'middle gray'

        pp.coordinates_to_pymol(cg)

    with tf.NamedTemporaryFile() as f:
        with tf.NamedTemporaryFile(suffix='.pml') as f1:
            f.write(pp.pymol_string())
            f.flush()

            pymol_cmd = 'hide all\n'
            pymol_cmd += 'run %s\n' % (f.name)
            pymol_cmd += 'show cartoon, all\n'
            pymol_cmd += 'bg white\n'

            f1.write(pymol_cmd)
            f1.flush()

            p = sp.Popen(['pymol', f1.name])
            out, err = p.communicate()

if __name__ == '__main__':
    main()

