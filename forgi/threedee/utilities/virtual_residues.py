from __future__ import division, unicode_literals

import numpy as np


def serialize_vres(vpos):
    """
    Output carthesian coordinates for each virtual residue, ordered by nucleotide number and seperated by space.

    The returned coordinates are relative to the element's coordinate system.

    :param vpos: A dictionary `{nucleotide_number: coordinates}`.
                 `nucleotide_number` is an integer with 0 referring to the first nucleotide of the element.
                 `coordinates` is a 3-element 1D numpy array
                 `nucleotide_numbers` should be consecutive integers.
    :returns: A string with space-seperated numbers.
              Groups of three numbers are x,y,z coordinates of one residue.
    """
    out_str = ""
    for k, v in sorted(vpos.items(), key=lambda x: x[0]):
        out_str += "{:.8f} {:.8f} {:.8f} ".format(v[0], v[1], v[2])
    return out_str


def parse_vres(parts):
    """
    Go back from a splitted serialized representation to a dictionary used for vres.

    :param parts: A list of strings which hold floats. The length has to be a multiple of 3.
                  Like a list gegerated by `serialize_vres(vpos).split()`
    :returns: A dictionary {nucleotide_number: coordinated}. See param `vpos` of the
              function `serialize_vres`
    """
    vpos = {}
    vbase = {}
    vsugar = {}
    vbackbone = {}
    to_fill = vpos
    start = 0
    for to_fill in [vpos, vbase, vsugar, vbackbone]:
        for i, j in enumerate(range(start, len(parts), 3)):
            if parts[j].startswith("v"):
                start = j + 1
                break
            to_fill[i] = np.array(list(map(float, parts[j:j + 3])))

    return vpos, vbase, vsugar, vbackbone
