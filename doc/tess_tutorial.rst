Tess Tutorial
==============
Introduction
------------
Tess is an extension of corgy capable of handling 3D data about RNA structures. It provides methods for extracting secondary structure as well as a coarse grain representation of 3D RNA structures.

Examples
--------

Extracting 2D structure from a 3D structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: 3d_to_2d.png
    :width: 800
    :align: center

To extract the base-pair information from a 3D structure stored in a PDB file, use the `pdb_to_ss_fasta.py` script::

    [pete@kat corgy]$ python examples/pdb_to_ss_fasta.py examples/1y26.pdb 
    >1y26
    CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG
    (((((((((...(((((((.......)))))))........((((((.......))))))..)))))))))


