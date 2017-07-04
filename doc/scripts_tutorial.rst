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
"coarse-grain file", "bpseq" and "pdb", use the corresponding scripts.

View simulated electron microscopy images of a cg-file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The script :code:`investigate_projection.py` opens a graphical window that lets you generate 
projected images at different projection angles and different resolution.

Get a representation of the coarse grain element names that correspond to dots and brackets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use::

    dotbracket_to_element_string.py -s -e 1y26_ss.dotbracket

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

    ../venv/bin/python graph_to_neato.py 1y26.cg | neato -Tpng -o /tmp/test.png 


