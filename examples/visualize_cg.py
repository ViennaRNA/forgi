#!/usr/bin/python

import sys
import subprocess as sp
import tempfile as tf

import forgi.threedee.model.coarse_grain as cmg
import forgi.utilities.debug as cud
import forgi.threedee.utilities.pdb as cup
import forgi.threedee.visual.pymol as cvp

from optparse import OptionParser

def main():
    usage = """
    ./visualize_cg.py cg_file

    Display the coarse-grain representation of a structure in pymol.
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')
    parser.add_option('-l', '--loops', dest='loops', default=True, action='store_false', help="Don't display the coarse-grain hairpin loops")

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    cg = cmg.from_file(args[0])
    pp = cvp.PymolPrinter()

    pp.add_loops = options.loops
    #cud.pv('cg.to_cg_string()')
    #sys.exit(1)
    pp.coordinates_to_pymol(cg)
    pp.print_text = False
    #pp.print_text = False
    #pp.output_pymol_file()

    with tf.NamedTemporaryFile() as f:
        with tf.NamedTemporaryFile(suffix='.pml') as f1:
            f.write(pp.pymol_string())
            f.flush()

            pymol_cmd = 'hide all\n'
            pymol_cmd += 'run %s\n' % (f.name)
            pymol_cmd += 'show cartoon, all\n'

            f1.write(pymol_cmd)
            f1.flush()

            p = sp.Popen(['pymol', f1.name])
            out, err = p.communicate()

if __name__ == '__main__':
    main()

