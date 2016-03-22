from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
import forgi.threedee.model.coarse_grain as ftmc
__metaclass__ = type #New style classes in Python 2.x

import forgi.projection.hausdorff as fph
import sys
import numpy as np
import scipy.misc
import scipy.ndimage
import matplotlib.pyplot as plt
import time
"""
A preliminary test before implementing the HausdorffEnergy in ernwin.
"""

class HausdorffEnergy(object):
    def __init__(self, img, scale):
        """
        :param img: A boolean, square 2D numpy.array
        :param scale: Int or float. How many Angstrom the side length of the image is.
        """
        super(HausdorffEnergy, self).__init__()
        self.ref_img=img
        self.ref_half=(scipy.ndimage.zoom(img, 0.5)>0.3)
        self.ref_quarter=(scipy.ndimage.zoom(img, 0.25)>0.3)
        self.scale=scale
        self.last_dir=None
    def shortname(self):
        return "HDF"
    def get_name(self):
        return "Hausdorff-Energy"
    def eval_energy(self, cg, background=None, nodes=None):
        start=time.time()
        s, i, self.last_dir, last_rot = fph.globally_minimal_distance(self.ref_quarter, 
                                                                      self.scale, cg, 
                                                                      virtual_atoms=False)
        end=time.time()
        td=end-start
        print("Global Opt took {}min {} sec".format(td//60, td%60))
        fig, ax=plt.subplots(2)
        ax[0].imshow(self.ref_quarter, interpolation="none", cmap='gray')
        ax[1].imshow(i, interpolation="none", cmap='gray')
        ax[0].set_title("Reference Global Optimization")
        ax[1].set_title("{} distance".format(s))
        plt.show()
        """ #The intermediate step with half resolution is not necessary!
        start=time.time()
        s, i, self.last_dir, last_rot = fph.locally_minimal_distance_one(self.ref_half, self.scale,
                                                                         cg, last_rot, self.last_dir, 100)
        end=time.time()
        td=end-start
        print("Local Opt took {}min {} sec".format(td//60, td%60))
        fig, ax=plt.subplots(2)
        ax[0].imshow(self.ref_half, interpolation="none", cmap='gray')
        ax[1].imshow(i, interpolation="none", cmap='gray')
        ax[0].set_title("Local Opt Half Reference")
        ax[1].set_title("{} distance".format(s))
        plt.show()"""
        start=time.time()
        score, img, self.last_dir, last_rot = fph.locally_minimal_distance_one(self.ref_img, self.scale, 
                                                                               cg, last_rot, self.last_dir, 100)
        end=time.time()
        td=end-start
        print("Local Opt took {}min {} sec".format(td//60, td%60))
        fig, ax=plt.subplots(2)
        ax[0].imshow(self.ref_img, interpolation="none", cmap='gray')
        ax[1].imshow(img, interpolation="none", cmap='gray')
        ax[0].set_title("Final Opt Reference")
        ax[1].set_title("{} distance".format(score))
        plt.show()
        return 10*score**2


if __name__=="__main__":
    cg=ftmc.CoarseGrainRNA(sys.argv[1])
    ref_img=scipy.misc.imread(sys.argv[2], flatten=True)
    energy=HausdorffEnergy(ref_img, 450)
    print(energy.eval_energy(cg))
