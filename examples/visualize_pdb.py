#!/usr/bin/python

from __future__ import print_function
import sys
import os.path as op
import subprocess as sp
import tempfile as tf
import logging
import warnings

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.visual.pymol as ftvp

from optparse import OptionParser

import numpy as np

def get_residue(chains, label):
    chain, number = label.split(":")
    return chains[chain][(' ', int(number), ' ')]

def main():
    logging.basicConfig()
    usage = """
    ./visualize_pdb.py pdb_file

    Display a pdb file along with its coarse grain representation in pymol.
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    parser.add_option('-s', '--secondary-structure',
                      dest='secondary_structure', default='',
                      help="Enter a dot-bracket string for the \
                      secondary structure of this model", type='str')
    parser.add_option('-x', '--text', dest='text', default=False, action='store_true', help="Add labels to the figure.")
    parser.add_option('-r', '--longrange', dest='longrange', default=False, action='store_true', help="Display long-range interactions")
    parser.add_option('-p', '--pseudoknots', dest='pseudoknots', default=False, action='store_true', help='Allow pseudoknots in the CG structure')
    #parser.add_option('', '--batch', dest='batch', default=False, action='store_true', help='Start pymol in batch mode') #Crashes currently
    parser.add_option('', '--virtual-residues', dest='virtual_residues', default=False, action='store_true', help='Display the virtual residues as spheres')
    parser.add_option('', '--color-residues', dest='color_residues', default=False, action='store_true', help="Color the residues according to the element they're in")

    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    parser.add_option('-l', '--loops', dest='loops', default=True, action='store_false', help="Don't display the coarse-grain hairpin loops")
    parser.add_option( '', '--chain', dest='chain', type='str', help="DWhat chain you like to display.")

    parser.add_option('-d', '--distance', dest='distance', default=None, help="Draw the lines between specified virtual residues")
    parser.add_option('-o', '--output', dest='output', default=None, help="Create a picture of the scene and exit", type='str')
    parser.add_option('', '--only-elements', dest='only_elements', default=None,
                     help='Display only these elements, separated by commas')
    parser.add_option('', '--thin-cylinders', dest='thin_cylinders', default=None,
                        help='Make coarse_grain RNA thinner')

    parser.add_option('-v', '--virtual-atoms', dest='virtual_atoms', default=False, action='store_true',
                      help='Display the virtual atoms')
    parser.add_option('', '--show-bp-detection', dest='show_bp_detection', type="str",
                      help="A comma-seperated list of residue numbers, "
                           "for which the forgi-basepair-detection algorithm will be illustrated.")
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    if not op.exists(args[0]):
        print("File doesn't exist: %s" % (args[0]), file=sys.stderr)
        sys.exit(1)


    if options.chain:
        chain_id=options.chain
    else:
        chain_id=None
    cgs = ftmc.CoarseGrainRNA.from_pdb(args[0], secondary_structure = options.secondary_structure.strip("\"'"),
                      remove_pseudoknots=not options.pseudoknots, load_chains=chain_id)
    cg = cgs[0]
    if len(cgs)>1:
        warnings.warn("The pdb has multiple connected components: {}. Showing the first".format([c.name for c in cgs]))
    pp = ftvp.PymolPrinter()
    if options.thin_cylinders:
        pp.cylinder_width=0.5
        pp.show_twists=False
    pp.display_virtual_residues = options.virtual_residues
    pp.virtual_atoms = options.virtual_atoms

    if options.only_elements is not None:
        orig_only_elements = options.only_elements.split(',')
        only_elements = set(orig_only_elements[::])
        for c in orig_only_elements:
            for e in cg.edges[c]:
                only_elements.add(e)
        pp.only_elements = only_elements

    pp.add_loops = options.loops
    pp.add_longrange = options.longrange
    #sys.exit(1)
    if not options.show_bp_detection:
        pp.coordinates_to_pymol(cg)
    pp.print_text = options.text
    print("virtual_residues:", options.virtual_residues, file=sys.stderr)
    #pp.print_text = False
    #pp.output_pymol_file()

    if options.only_elements is not None:
        pp.only_elements = options.only_elements.split(',')

    with tf.NamedTemporaryFile(mode="w+") as f:
        with tf.NamedTemporaryFile(suffix='.pml',mode="w+") as f1:
            with tf.NamedTemporaryFile(suffix='.pdb',mode="w+") as f2:
                chains = cg.chains
                #chain = ftup.renumber_chain(chain)
                ftup.output_multiple_chains(chains.values(), f2.name)
                f2.flush()

                # display the distances between nucleotides
                if options.distance is not None:
                    for dist_pair in options.distance.split('*'):
                        fr,to = dist_pair.split(',')
                        try:
                            vec1 = get_residue(chains, fr)["C1'"].get_vector().get_array()
                            vec2 = get_residue(chains, to)["C1'"].get_vector().get_array()
                        except KeyError:
                            # Rosetta produces atoms with non-standard names
                            vec1 = get_residue(chains, fr)["C1*"].get_vector().get_array()
                            vec2 = get_residue(chains, to)["C1*"].get_vector().get_array()

                        pp.add_dashed(vec1, vec2, width=1.2)


                if options.show_bp_detection:
                    r_key1, r_key2 = options.show_bp_detection.split(",")
                    res1 = get_residue(chains,r_key1)
                    res2 = get_residue(chains, r_key2)
                    points = ftup._get_points(res1, res2)
                    for pair in points[1:]:
                        if ftuv.vec_distance(pair[0], pair[1])<ftup.HBOND_CUTOFF:
                            c="green"
                        else:
                            c="purple"
                        pp.add_dashed(pair[0], pair[1], width=0.2, color=c)
                    flattened_points = [ p for pair in points for p in pair]
                    coplanar_indices, ctr, normal_vec = ftup._coplanar_point_indices(*flattened_points)
                    for i,p in enumerate(flattened_points):
                        if i in coplanar_indices:
                            c = "green"
                        else:
                            c="purple"
                        pp.add_sphere(p, color=c, width=0.4, text="")

                    # Plane: n*P=n*A or n*(P-A)=0

                    # The two reference points not involved in base-pairing
                    d = np.dot(normal_vec, ctr)
                    x1 = points[0][0][0]
                    y1 = points[0][0][1]
                    # ax+by+cz=d
                    z1 = (d-(normal_vec[0]*x1+normal_vec[1]*y1))/normal_vec[2]
                    x3 = points[0][1][0]
                    y3 = points[0][1][1]
                    z3 = (d-(normal_vec[0]*x3+normal_vec[1]*y3))/normal_vec[2]
                    FACTOR=1.2
                    # Get the normal to the connection of the two points in the plane
                    conn = np.array([ x1-x3, y1-y3, z1-z3 ])
                    vect = np.cross(conn, normal_vec)
                    p1 = np.array([x1,y1,z1])-FACTOR*conn
                    p3 = np.array([x3,y3,z3])+FACTOR*conn
                    l = ftuv.vec_distance(p1,p3)/2
                    vect/=ftuv.magnitude(vect)
                    vect*=l
                    corners = [ p1, ctr+vect, p3, p3, ctr-vect, p1]
                    pp.boxes.append((corners, 'yellow'))
                print(pp.pymol_string(), file=f)
                f.flush()

                pymol_cmd = 'hide all\n'
                pymol_cmd += 'show cartoon, all\n'
                pymol_cmd += 'set cartoon_ring_mode\n'
                pymol_cmd += 'set cartoon_tube_radius, .3\n'

                if options.only_elements is not None:
                    pymol_cmd += "hide all\n"

                    for constraint in only_elements:
                        color = pp.get_element_color(constraint)

                        for r in cg.define_residue_num_iterator(constraint, seq_ids=True):
                            pymol_cmd += "show sticks, resi %r\n" % (r[1])
                            pymol_cmd += "color %s, resi %r\n" % (color, r[1])

                if options.color_residues:
                    for d in cg.defines:
                        color = pp.get_element_color(d)

                        for r in cg.define_residue_num_iterator(d, seq_ids=True):
                            pymol_cmd += "color %s, resi %r\n" % (color, r[1])


                pymol_cmd += 'run %s\n' % (f.name)
                pymol_cmd += 'bg white\n'
                pymol_cmd += 'clip slab, 10000\n'
                pymol_cmd += 'orient\n'

                if options.output is not None:
                    pymol_cmd += 'ray\n'
                    pymol_cmd += 'png %s\n' % (options.output)
                    #pymol_cmd += 'quit\n' #This would lead to an error message

                f1.write(pymol_cmd)
                f1.flush()

                #if options.batch:
                #    p = sp.Popen(['pymol', '-cq', f2.name, f1.name], stdout=sp.PIPE, stderr=sp.PIPE)
                if True: #else:
                    p = sp.Popen(['pymol', f2.name, f1.name], stdout=sp.PIPE, stderr=sp.PIPE)

                out, err = p.communicate()

if __name__ == '__main__':
    main()
