.. _installation:

Installation
============

Download
--------

The source code for the forgi module can be found at github:

https://github.com/pkerpedjiev/forgi

To install, just run the usual::

    sudo python setup.py install

Or install to a particular directory::

    python setup.py install --prefix=/home/pete/local

Dependencies
------------

Some standard python packages are required:

* numpy + scipy
* the future package
* networkx


Additional dependencies
~~~~~~~~~~~~~~~~~~~~~~~

The library can be installed with some of the additional dependencies missing, but will 
raise an error, if functions with missing dependencies are called.

These additional dependencies include:

* matplotlib (Needed for `forgi.threedee.visual`, `forgi.visual` and `forgi.projection`)
* Biopython (Needed for `forgi.threedee`)
* pylab (Needed for `forgi.visual.mplotlib`, )

MC-Annotate
~~~~~~~~~~~

To convert PDB files into CoarseGrainRNA objects, MC-Annotate is required.
It is available at http://major.iric.ca/MajorLabEn/MC-Tools.html
