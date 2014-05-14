RNA 3D Structure Using forgi.threedee
==============
Introduction
------------
`forgi.threedee` is an extension of forgi capable of handling 3D data about RNA structures. It provides methods for extracting secondary structure as well as creating a coarse grain representation of 3D RNA structures.

Requirements
------------

MC_Annotate_ is required for annotating the nucleotide interactions.

.. _MC_Annotate: http://www.major.iric.ca/MajorLabEn/MC-Tools.html

Examples
--------

Extracting 2D structure from a 3D structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: 3d_to_2d.png
    :width: 800
    :align: center

To extract the base-pair information from a 3D structure stored in a PDB file, use the `pdb_to_ss_fasta.py` script. By default, this script operates on the largest RNA chain in the provided pdb file (although a specific chain can be specified using the `-c` option).::

    [pete@kat forgi]$ python examples/pdb_to_ss_fasta.py examples/1y26.pdb 
    >1y26
    CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG
    (((((((((...(((((((.......)))))))........((((((.......))))))..)))))))))

The reported 2D structure is extracted from the annotations of MC-Annotate [1]. Pseudoknots are removed using the knotted2nested.py script which is included with the generous permition of Dr. Sandra Smit [2]. Users making use of the 3D-to-2D facilities of the ``forgi.threedee`` package should cite the articles listed in the citations section at the bottom of this page.

Creating a Coarse Grain 3D Representation of an RNA Molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can imagine an RNA molecule as a collection of helices and ... not helices (hairpins, interior loops, multiloops, etc..) as described in the :ref:`forgi_graph_tutorial`. By creating an abstraction for the canonical helixes and representing them as cylinders, we can create a coarse grain representaiton of the RNA molecule::

    [pkerp@plastilin forgi]$ python examples/pdb_to_cg.py ~/projects/ernwin/fess/structures/2l94.pdb 
    name 2l94
    length 45
    seq GGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAGGGAAUCUUCCC
    define f1 0 1
    define i0 34 38 8 9
    define s1 9 19 24 34
    define s0 1 8 38 45
    define t1 45 46
    define h0 19 24
    connect s1 i0 h0
    connect s0 f1 i0 t1
    coord i0 -3.16011516235 3.53600076719 3.59573158587 -10.6370223437 3.13588043498 -0.447155058784
    coord s1 -3.16011516235 3.53600076719 3.59573158587 21.5002000702 -2.08161219996 0.501285407577
    coord s0 -25.7014473529 -5.72949740183 -4.01197045617 -10.6370223437 3.13588043498 -0.447155058784
    coord h0 21.5002000702 -2.08161219996 0.501285407577 23.5599994659 6.76999998093 -5.87900018692
    twist s1 -0.00919853915562 0.451186349593 -0.892382353489 -0.153976361554 -0.136900113113 -0.978544653612
    twist s0 -0.532923452202 0.816983142634 0.220297840992 0.430691441325 -0.408188889295 -0.804914102887

The lines beginning with `coord` indicate the positions of the helices. They are calculated by fitting a cylinder to each Watson-crick paired region of the 3D structure.

Visualizing the Coarse Grain Helical Representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The coarse grain representation shown in the prior example can be visualized using Pymol_ and the ``visualize_pdb.py`` script. The green cylinder below represent canonical helical regions, red represents multiloops and blue - hairpin loops::

    python examples/visualize_pdb.py ~/projects/ernwin/fess/structures/1y26.pdb

.. image:: 1y26_visualized.png
    :width: 400
    :align: center

.. _Pymol: http://www.pymol.org/

If you just have the coarse-grain file, then use the ``visualize_cg.py`` script::

    python examples/visualize_cg.py examples/1y26.cg

.. image:: 1y26_coarse.png
    :width: 400
    :align: center

Get A Description of a Coarse-Grain Stem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `get_stem_stats` function returns a `forgi.threedee.model.StemStat` structure which contains some information about a particular stem, such as how many base pairs it has, how long it is (in Angstroms) and how much its helix twists from the start to the end. It also stores information about which nucleotides it is composed of (its `define`). 

Using the structure 2MIS as an example::

    >>> import forgi.threedee.model.coarse_grain as ftmc
    >>> cg = ftmc.from_pdb('test/forgi/threedee/data/2mis.pdb', intermediate_file_dir='tmp')
    >>> print cg.get_stem_stats('s0')
    pdb_name: 2mis_A bp_length: 6 phys_length: 14.735000 twist_angle: 2.822735 define: 1 6 21 26

This indicates that the first stem in the structure ('s0'), composed of the nucleotides 1 - 6 and 21 - 26 has a length of 14.735 Angstroms and a twist of 2.82 radians. It is composed of 6 base pairs and comes from a structure named `2mis_A`.

Citations
~~~~~~~~~

[1] *Gendron P, Lemieux S, Major F(2001)*. **Quantitative analysis of nucleic acid three-dimensional structures.** J Mol Biol 308:919–936.

[2] *Sandra Smit, Kristian Rother, Jaap Heringa, and Rob Knight*.
**From knotted to nested RNA structures: a variety of computational methods for pseudoknot removal.**
RNA (2008) 14(3):410-416.

