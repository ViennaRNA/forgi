# pylint: disable=W0612

"""
forgi - an RNA structure library
================================

Forgi is a python library for analysing and manipulating RNA secondary and
tertiary structures based on secondary structure elements.

The libraries core object, which is returned by the load_rna function,
is the `forgi.graph.bulge_graph.BulgeGraph` object (2D structure) and its
subclass the `forgi.threedee.model.coarse_grain.CoarseGrainRNA` object.

Citing
------

If you use forgi for a scientific publication, please cite this library
as specified in the CREDITS section.
"""

__author__ = "Bernhard C. Thiel, Peter Kerpedjiev"
__copyright__ = "Copyright 2012 - 2019"
__license__ = "GNU GPL v 3.0"
__version__ = "2.2.3"
__maintainer__ = "University of Vienna, TBI"


import os

from forgi.utilities.commandline_utils import load_rna


def data_file(fname):
    """Return the path to a data file of ours."""
    return os.path.join(os.path.split(__file__)[0], fname)
