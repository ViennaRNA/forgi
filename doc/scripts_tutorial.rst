.. _forgi_scripts_tutorial:

Useful scripts coming with forgi
================================

Forgi comes with a few scripts in the example folder.
Here some of the most useful are documented.

Visualization scripts
~~~~~~~~~~~~~~~~~~~~~

There are scripts for visualizing PDB and coarse-grain RNA structures
(See :doc:`graph_tutorial` and :doc:`threedee_tutorial`
for more information on Coarse Grain Files).

To view a pdb file use::

    visualize_pdb.py 1jj2.pdb

To view a cg file use::

    visualize_cg.py 1y26.cg

Use :code:`visualize_cg.py -x 1y26.cg` to show coarse-grained element names.

Use :code:`visualize_cg.py --virtual-atoms 1y26.cg` to show virtual atoms for the backbone.

Use :code:`visualize_cg.py --sidechain-atoms 1y26.cg` to show all virtual atoms.

File format conversions
~~~~~~~~~~~~~~~~~~~~~~~

To convert files between the file formats "fasta with secondary structure", "dotbracket string",
"coarse-grain file", "bpseq" and "pdb", use the :code:`rnaConvert.py` script::

    rnaConvert.py 1jj2.pdb -T forgi

The :code:`-T` option specifies the output file format.

Use :code:`-T forgi` to create a "*.bg" or "*.cg" file, :code:`-T fasta` for a
fasta-file with secondary structure, :code:`-T bpseq` for a bpseq-file and
:code:`-T dotbracket` to only output a dotbracket string.

You can use `--filename OUTNAME` to write to a file instead of STDOUT. Note that the file extension
will be appended automatically.
If the input-file contains multiple RNA molecules (connected components),
a file will be created for each of them.


View simulated electron microscopy images of a cg-file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The script :code:`investigate_projection.py` opens a graphical window that lets you generate
projected images at different projection angles and different resolution.

Get a representation of the coarse grain element names that correspond to dots and brackets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use::

    rnaConvert.py 1y26_ss.dotbracket -T element_string

This will print::

    (((((((((...(((((((.......)))))))........((((((.......))))))..)))))))))
    sssssssssmmmssssssshhhhhhhsssssssmmmmmmmmsssssshhhhhhhssssssmmsssssssss
    00000000000011111110000000111111122222222222222111111122222211000000000

The numbers indicate the numbers of the coarse grained elements. As expected, the first column
holds the values '(','s','0' indicating that this opening bracket belongs to the stem "s0".
This is especially useful when writing tests for the forgi library or code depending on it.

Render the graph-representation of a Bulge Graph as an image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the following command, which requires neato, to generate a png image
of the graph represenation of the RNA::

    rnaConvert.py 1y26.cg -T neato8 | neato -Tpng -o /tmp/test.png


Update cg-files created with earlier versions of forgi
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you use CoarseGrainRNA files (*.cg/*.coord files) that were created before release of version 1.0,
you might get an :code:`AssertionError` saying that the twists are inconsistent. Due to a bug in
earlier versions of forgi, some twist vectors could sometimes deviate by a few degrees from
being orthogonal to the stem vector. To use such broken *.cg files with newer versions of forgi,
you can use the example script: :code:`fix_twists.py`::

A list of all scripts
=====================

*  :code:`all_atom_positions.py` and :code:`average_atom_positions.py` were used to generate the file
`forgi/threedee/data/average_atom_positions.json`, which comes with forgi.

*  `burial.py`: Characterize a pdb file.

* `compare_RNA.py`: Calculate RMSD and ACC between 2 3D structures.
