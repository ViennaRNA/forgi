#!/usr/bin/python

import sys
import os.path as op
import subprocess as sp
import tempfile as tf

import forgi.threedee.model.coarse_grain as ftmc
import forgi.utilities.debug as fud
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.visual.pymol as ftvp

from optparse import OptionParser

def main():
    usage = """
    ./visualize_pdb.py pdb_file

    Display a pdb file along with its coarse grain representation in pymol.
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    parser.add_option('-s', '--secondary-structure', 
                      dest='secondary_structure', default='', 
                      help="Enter a dot-bracket string for the \
                      secondary structure of this model", type=str)
    parser.add_option('-x', '--text', dest='text', default=False, action='store_true', help="Add labels to the figure.")
    parser.add_option('-r', '--longrange', dest='longrange', default=False, action='store_true', help="Display long-range interactions")
    parser.add_option('-c', '--constraints', dest='constraints', default=None, help="Only visualize the elements passed as parameters", type='str')
    parser.add_option('-p', '--pseudoknots', dest='pseudoknots', default=False, action='store_true', help='Allow pseudoknots in the CG structure')

    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    parser.add_option('-l', '--loops', dest='loops', default=True, action='store_false', help="Don't display the coarse-grain hairpin loops")

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    if not op.exists(args[0]):
        print >>sys.stderr, "File doesn't exist: %s" % (args[0])
        sys.exit(1)

    cg = ftmc.from_pdb(args[0], options.secondary_structure.strip("\"'"),
                      remove_pseudoknots=not options.pseudoknots)
    pp = ftvp.PymolPrinter()

    if options.constraints is not None:
        orig_constraints = options.constraints.split(',')
        constraints = set(orig_constraints[::])
        for c in orig_constraints:
            for e in cg.edges[c]:
                constraints.add(e)
        pp.constraints = constraints

    pp.add_loops = options.loops
    pp.add_longrange = options.longrange
    #sys.exit(1)
    pp.coordinates_to_pymol(cg)
    pp.print_text = options.text
    #pp.print_text = False
    #pp.output_pymol_file()

    with tf.NamedTemporaryFile() as f:
        with tf.NamedTemporaryFile(suffix='.pml') as f1:
            with tf.NamedTemporaryFile(suffix='.pdb') as f2:
                # extract just the biggest chain and renumber it so
                # the nucleotides start at 1
                chain = ftup.get_biggest_chain(args[0])
                #chain = ftup.renumber_chain(chain)
                ftup.output_chain(chain, f2.name)
                f2.flush()

                f.write(pp.pymol_string())
                f.flush()

                pymol_cmd = 'hide all\n'
                pymol_cmd += 'show cartoon, all\n'
                pymol_cmd += 'set cartoon_ring_mode\n'
                pymol_cmd += 'set cartoon_tube_radius, .3\n'

                if options.constraints is not None:
                    pymol_cmd += "hide all\n"

                    for constraint in constraints:
                        color = pp.get_element_color(constraint)

                        for r in cg.define_residue_num_iterator(constraint):
                            pymol_cmd += "show sticks, resi %r\n" % (r)
                            pymol_cmd += "color %s, resi %r\n" % (color, r)

                pymol_cmd += 'run %s\n' % (f.name)

                f1.write(pymol_cmd)
                f1.flush()

                p = sp.Popen(['pymol', f2.name, f1.name])
                out, err = p.communicate()

if __name__ == '__main__':
    main()

