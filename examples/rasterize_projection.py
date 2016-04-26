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

WIDTH=350
if __name__=="__main__":
    for filename in sys.argv[1:]:
        cg=ftmc.CoarseGrainRNA(filename)
        try:
            proj=fpp.Projection2D(cg, project_virtual_atoms=True)
        except ValueError:
            a=random.random()
            b=random.random()
            c=random.random()
            proj=fpp.Projection2D(cg, project_virtual_atoms=True, proj_direction=[a,b,c])
        rot=random.randrange(0,3600)/10
        proj.rotate(rot)
        box=fph.get_box(proj, WIDTH)
        img, _= proj.rasterize(700, box)
        scipy.misc.imsave(filename+".rot{}.width{}.png".format(rot, WIDTH), img)
