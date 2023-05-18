Forgi Library for RNA Secondary Structure Analysis
==================================================

.. image:: https://travis-ci.org/ViennaRNA/forgi.svg?branch=master
    :target: https://travis-ci.org/ViennaRNA/forgi



.. image:: https://zenodo.org/badge/6291657.svg
   :target: https://zenodo.org/badge/latestdoi/6291657

.. image:: https://anaconda.org/bioconda/forgi/badges/version.svg   
   :target: https://anaconda.org/bioconda/forgi

Full documentation: https://viennarna.github.io/forgi

Citing
======

If your use of forgi results in an academic publication, please cite:

Thiel BC, Beckmann IK, Kerpedjiev P and Hofacker IL. 3D based on 2D: Calculating helix angles and stacking patterns using forgi 2.0, an RNA Python library centered on secondary structure elements. [version 2; peer review: 2 approved]. *F1000Research* 2019, **8**:287
(https://doi.org/10.12688/f1000research.18458.2) 

For the pseudoknot-removal code and external PDB annotation tools, see the citations therein.

Installation
============

:code:`pip install forgi` will work in many cases. However, as forgi contains a small part of compiled code, this can fail,
especially on operating systems/ python versions for which no manylinux wheel is available.

If compilation fails, you should first try to install Cython, python development headers and a compiler. E.g. on Ubuntu
:code:`sudo apt install gcc g++ python3-dev`. On other distributions/ when targeting non-default python versions, the package cpoyuld also be called python3-devel or python3.11-dev or something similar.

Locally build the documentation 
===============================

:code:`sphinx-build -b html doc OUTPUTFOLDER`

