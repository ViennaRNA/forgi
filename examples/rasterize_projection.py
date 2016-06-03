from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

__metaclass__ = type #New style classes in Python 2.x

import forgi.threedee.model.coarse_grain as ftmc
import forgi.projection.projection2d as fpp
import forgi.projection.hausdorff as fph

import sys, random
import numpy as np
import scipy.misc
import scipy.ndimage
import argparse

WIDTH=350

def get_parser():
    """
    Here all commandline & help-messages arguments are defined.

    :returns: an instance of argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    #Argument(s)
    parser.add_argument('cgfiles', nargs='+', help='One or more *.cg/*.coord files holding the RNA to plot.')
    parser.add_argument('-w', '--width', type=int ,help='Width in Angstrom', default = 300)
    parser.add_argument('-d', '--dpi', type=int, help='Number of pixels in the image', default = 50)
    parser.add_argument('-r', '--rotate', type=int, help='Degrees to rotate in plane', default=99999999999)
    parser.add_argument('-n', type=int, help='Number of images to produce.', default = 1)
    return parser

parser = get_parser()
if __name__=="__main__":
    args = parser.parse_args()
    for filename in args.cgfiles:
        for n in range(args.n):
            cg=ftmc.CoarseGrainRNA(filename)
            try:
                proj=fpp.Projection2D(cg, project_virtual_atoms=True)
            except ValueError:
                a=random.random()
                b=random.random()
                c=random.random()
                proj=fpp.Projection2D(cg, project_virtual_atoms=True, proj_direction=[a,b,c])
            if args.rotate==99999999999:
                rot=random.randrange(0,3600)/10
            else:
                rot=args.rotate
            proj.rotate(rot)
            box=fph.get_box(proj, args.width)
            img, _= proj.rasterize(args.dpi, box)
            outname=filename+".dpi{}.width{}.{}.png".format(args.dpi, args.width, n)
            scipy.misc.imsave(outname, img, "png")
            print("File {} written".format(outname))
