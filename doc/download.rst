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

* numpy + scipy + pandas
* the future package
* networkx


Further more logging_exceptions is required, which we published separately on Pypi

    pip install logging_exceptions



Additional dependencies
~~~~~~~~~~~~~~~~~~~~~~~

The library can be installed with some of the additional dependencies missing, but will
raise an error, if functions with missing dependencies are called.

These additional dependencies include:

* matplotlib (Needed for `forgi.threedee.visual`, `forgi.visual` and `forgi.projection`)
* Biopython (Needed for `forgi.threedee`)
* pylab (Needed for `forgi.visual.mplotlib`, )
* beautifulsoap4 (Needed for fetching informations about modified residues in pdb files)
* appdirs (Needed for `forgi.threedee`)

Either MC-Annotate or DSSR
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To convert PDB files into CoarseGrainRNA objects, either MC-Annotate or DSSR is required.

MC-Annotate is available at http://major.iric.ca/MajorLabEn/MC-Tools.html

DSSR is available at http://home.x3dna.org/

You can create a configuration file to decide which of the two programs you want to use.
If no configuration file exists, you will be prompted to choose the program the first time you
convert a PDB file. Your choice will be recorded in a newly generated configuration file.

We use the appdirs package to define the system specific location of the
configuration file. Under Linux it would default to ~/.config/forgi.
