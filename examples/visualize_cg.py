#!/usr/bin/python

from __future__ import print_function
from __future__ import division
import numpy as np
import sys
import subprocess as sp
import tempfile as tf
import time

import forgi.threedee.model.coarse_grain as cmg
import forgi.utilities.debug as fud
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.pdb as cup
import forgi.threedee.model.similarity as ftur
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.visual.pymol as cvp

from optparse import OptionParser


def align_cgs(cgs):
    '''
    Align each coarse grain RNA to the first one.

    The points representing each coarse grain RNA molecule
    will be the virtual residues.
    
    @param cgs: A list of CoarseGrainRNA structures.
    @return: Nothing, the cgs are modified in place
    '''
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
    num_args = 1
    parser = OptionParser(usage=usage)

    # parser.add_option('-u', '--useless', dest='uselesss',
    # default=False, action='store_true', help='Another useless option')
    parser.add_option('-g', '--highlight', dest='highlight', default=None, help="Highlight some elements", type='str')
    parser.add_option('-o', '--output', dest='output', default=None, help="Create a picture of the scene and exit",
                      type='str')
    parser.add_option('-r', '--longrange', dest='longrange', default=False, action='store_true',
                      help="Display long-range interactions")
    parser.add_option('-l', '--loops', dest='loops', default=True, action='store_false',
                      help="Don't display the coarse-grain hairpin loops")
    parser.add_option('-c', '--cones', dest='cones', default=False, action='store_true',
                      help="Display cones that portrude from the stems")
    parser.add_option('-x', '--text', dest='text', default=False, action='store_true', help="Add labels to the figure.")
    parser.add_option('-a', '--align', dest='align', default=False, action='store_true',
                      help='Align all of the structures with the first')
    parser.add_option('-e', '--encompassing-stems', dest='encompassing_stems', default=False, action='store_true',
                      help='Show the big stems that encompass the colinear ones.')
    parser.add_option('-v', '--virtual-atoms', dest='virtual_atoms', default=False, action='store_true',
                      help='Display the virtual atoms')
    parser.add_option('-d', '--distance', dest='distance', default=None,
                      help="Draw the lines between specified virtual residues")
    parser.add_option('-t', '--residue-distance', dest='residue_distance', default=None,
                      help="Draw a line between residue distances")
    parser.add_option('-b', '--basis', dest='basis', default=False, action='store_true',
                      help='Display the coordinate basis of each element')
    parser.add_option('', '--stem-color', dest='stem_color', default='green',
                      help='The default color in coarse-grain drawings')
    parser.add_option('', '--multiloop-color', dest='multiloop_color', default='red',
                      help='The default color in coarse-grain drawings')
    parser.add_option('', '--batch', dest='batch', default=False, action='store_true', help='Start pymol in batch mode')
    parser.add_option('', '--sidechain-atoms', dest='sidechain_atoms', default=False, action='store_true',
                      help='Include the sidechain atoms. Automatically enables --virtual-atoms')
    parser.add_option('', '--rainbow', dest='rainbow', default=False, action='store_true',
                      help='Color each of the nucleotide positions (i.e. average atoms) according to the colors of \
                      the rainbow and their position')
    parser.add_option('', '--only-elements', dest='only_elements', default=None, help='Display only these elements '
                                                                                      'element names should be '
                                                                                      'separated by commas')
    parser.add_option('', '--color-gradual', dest='color_gradual', default=None, help='Color the specified elements'
                                                                                      'gradually from one to the other, example (i1,i4,m1)', type='str')

    (options, args) = parser.parse_args()

    print("hi")
    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)
    print("hi1")

    pp = cvp.PymolPrinter()
    pp.stem_color = options.stem_color
    pp.multiloop_color = options.multiloop_color
    pp.add_loops = options.loops
    pp.draw_cones = options.cones
    # sys.exit(1)
    pp.add_longrange = options.longrange
    pp.print_text = options.text
    pp.encompassing_stems = options.encompassing_stems
    pp.virtual_atoms = options.virtual_atoms
    pp.sidechain_atoms = options.sidechain_atoms
    pp.basis = options.basis
    pp.rainbow = options.rainbow

    if options.only_elements is not None:
        pp.only_elements = options.only_elements.split(',')

    cgs = []
    for a in args:
        cgs += [cmg.CoarseGrainRNA(a)]

    if options.align:
        align_cgs(cgs)

    if options.color_gradual is not None:
        pp.element_specific_colors = dict()
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap('coolwarm')

        for d in cgs[0].defines:
            pp.element_specific_colors[d]= 'black'

        to_color_nodes = options.color_gradual.split(',')
        for i,node in enumerate(to_color_nodes):
            print(node, cmap(i / float(len(to_color_nodes))))
            pp.element_specific_colors[node] = cmap(i / float(len(to_color_nodes)))

    for i, cg in enumerate(cgs):
        if i > 0:
            pp.color_modifier = .3
            #pp.override_color = 'middle gray'

        pp.coordinates_to_pymol(cg)


    # highlight things in purple
    if options.highlight is not None:
        for s in options.highlight.split(','):
            fud.pv('s')
            pp.add_twists = False
            pp.add_stem_like(cg, s, color='purple', width=3.)

    # display the distances between nucleotides
    if options.distance is not None:
        virtual_atoms = ftug.virtual_atoms(cg, sidechain=False)

        for dist_pair in options.distance.split(':'):
            fud.pv('dist_pair')
            fr, to = dist_pair.split(',')

            fr = int(fr)
            to = int(to)

            pp.add_dashed(virtual_atoms[fr]["C1'"], virtual_atoms[to]["C1'"], width=1.2)

    if options.residue_distance is not None:
        dist_pair = options.residue_distance
        fr, to = dist_pair.split(',')

        fr = int(fr)
        to = int(to)

        node1 = cg.get_node_from_residue_num(to)
        node2 = cg.get_node_from_residue_num(fr)

        pos1, len1 = cg.get_position_in_element(to)
        pos2, len2 = cg.get_position_in_element(fr)

        #fud.pv('node1, node2, pos1, pos2')

        vec1 = cg.coords[node1][1] - cg.coords[node1][0]
        vec2 = cg.coords[node2][1] - cg.coords[node2][0]

        #mid1 = (cg.coords[node1][0] + cg.coords[node1][1]) / 2
        #mid2 = (cg.coords[node2][0] + cg.coords[node2][1]) / 2

        mid1 = cg.coords[node1][0] + pos1 * (vec1 / len1)
        mid2 = cg.coords[node2][0] + pos2 * (vec2 / len2)

        pp.add_sphere(mid1, 'green', width=2)
        pp.add_sphere(mid2, 'red', width=2)
        

    with tf.NamedTemporaryFile() as f:
        with tf.NamedTemporaryFile(suffix='.pml') as f1:
            f.write(pp.pymol_string())
            f.flush()

            pymol_cmd = 'hide all\n'
            pymol_cmd += 'run %s\n' % (f.name)
            pymol_cmd += 'show cartoon, all\n'
            pymol_cmd += 'bg white\n'
            pymol_cmd += 'clip slab, 10000\n'
            pymol_cmd += 'orient\n'

            if options.output is not None:
                pymol_cmd += 'ray\n'
                pymol_cmd += 'png %s\n' % (options.output)
                pymol_cmd += 'quit\n'

            f1.write(pymol_cmd)
            f1.flush()

            print("f1.name:", f1.name)

            if options.batch:
                p = sp.Popen(['pymol', '-cq', f1.name], stdout=sp.PIPE, stderr=sp.PIPE)
            else:
                p = sp.Popen(['pymol', f1.name], stdout=sp.PIPE, stderr=sp.PIPE)

            out, err = p.communicate()
            print("err:", err, file=sys.stderr)

if __name__ == '__main__':
    main()

