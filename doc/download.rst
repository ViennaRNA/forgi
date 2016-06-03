Download
========

The source code for the forgi module can be found at github:

https://github.com/pkerpedjiev/forgi

To install, just run the usual::

    sudo python setup.py install

Or install to a particular directory::

    python setup.py install --prefix=/home/pete/local

Dependencies
============

The library can be installed with some of the dependencies missing, but will 
raise an error, if functions with missing dependencies are called.

Some standard python packages are required:

* numpy + scipy
* matplotlib

MC-Annotate
-----------

To convert PDB files into CoarseGrainRNA objects, MC-Annotate is required.
It is available at http://major.iric.ca/MajorLabEn/MC-Tools.html
