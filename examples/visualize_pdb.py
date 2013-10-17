#!/usr/bin/python

import sys
import subprocess as sp
import tempfile as tf

import corgy.model.coarse_grain as cmg
import corgy.utilities.debug as cud
import corgy.visual.pymol as cvp

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
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    cg = cmg.from_pdb(args[0], options.secondary_structure.strip("\"'"))

    pp = cvp.PymolPrinter()
    pp.coordinates_to_pymol(cg)
    pp.print_text = False
    #pp.print_text = False
    #pp.output_pymol_file()

    with tf.NamedTemporaryFile() as f:
        with tf.NamedTemporaryFile(suffix='.pml') as f1:
            cud.pv('pp.pymol_string()')

            f.write(pp.pymol_string())
            f.flush()

            pymol_cmd = 'hide all\n'
            pymol_cmd += 'run %s\n' % (f.name)
            pymol_cmd += 'show cartoon, all\n'

            f1.write(pymol_cmd)
            f1.flush()

            p = sp.Popen(['pymol', args[0], f1.name])
            out, err = p.communicate()

if __name__ == '__main__':
    main()

