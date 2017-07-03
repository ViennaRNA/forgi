#!/usr/bin/python

from __future__ import division
from builtins import range

import numpy as np
import os.path as op
import subprocess as sp
import sys
import tempfile as tf
from optparse import OptionParser

import forgi.threedee.utilities.average_atom_positions as ftua
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.visual.pymol as ftvp

import forgi.utilities.debug as fud

def main():
    usage = """
    python plot_average_atom_positions.py

    Currently works for certain multiloops (hardcoded).
    """

    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')
    parser.add_option('-a', '--all', dest='all_entries_file', default=None, help='Use a file containing the positions of each entry for this cg element.', type='str')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    pp = ftvp.PymolPrinter()

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap('gist_rainbow')
    stem_len = 3

    for aname in ftup.nonsidechain_atoms:
        conn_types = [2]
        for i in range(3):
            all_coords = []
            for conn_type in conn_types:
                elem_type = 'm'
                dimensions = '3 1000'
                identifier = "%s %s %d %d %s" % (elem_type,
                                              dimensions,
                                              conn_type,
                                              i, aname)

                if options.all_entries_file is not None:
                    import imp
                    fn = op.expanduser(options.all_entries_file)
                    #fn1 = op.expanduser(op.splitext(options.all_entries_file)[0])
                    xftua = imp.load_source('xftua', fn)
                    coords = xftua.all_atom_poss[identifier]
                    for c in coords:
                        pp.add_sphere(c, width=0.2, color_rgb = cmap(i / float(stem_len)))

                    all_coords += coords
                    avg_coord = np.mean(coords, axis=0)
                    pp.add_sphere(avg_coord, width=0.6, color_rgb = cmap(i / float(stem_len)))
                else:
                    coords = ftua.avg_atom_poss[identifier]
                    pp.add_sphere(coords, width=0.2, color_rgb = cmap(i / float(stem_len)))

            avg_coord = np.mean(all_coords, axis=0)
            pp.add_sphere(avg_coord, width=0.9, color_rgb = cmap(i / float(stem_len)))

    # Add an axis
    pp.add_segment([0,0,0], [5,0,0], "purple", 0.5)

    with tf.NamedTemporaryFile() as f:
        with tf.NamedTemporaryFile(suffix='.pml') as f1:
            f.write(pp.pymol_string())
            f.flush()

            pymol_cmd = 'hide all\n'
            pymol_cmd += 'run %s\n' % (f.name)
            pymol_cmd += 'show cartoon, all\n'
            pymol_cmd += 'bg white\n'
            pymol_cmd += 'clip slab, 10000\n'

            f1.write(pymol_cmd)
            f1.flush()

            p = sp.Popen(['pymol', f1.name])
            out, err = p.communicate()

if __name__ == '__main__':
    main()

