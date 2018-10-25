.. _forgi_utilities_tutorial:

RNA Secondary Structure Utilities Using the forgi.utilities Submodule
=====================================================================

Introduction
------------

forgi is a library for manipulating RNA as a graph-like
structure. It provides classes for reading, storing, manipulating and
outputting RNA structures. It classifies and stores RNA secondary
structure in a data structure that quickly allows one to determine
which structural  elements are connected and which nucleotides they contain.

However, some taske ar ebetter performed on the dotbracket structure or
the pair-table representation of an RNA. For this, the `forgi.utilities` module
provides useful functions.

Dot-bracket to pair-table
-------------------------

Converting a secondary structure in dot-bracket format to a pairtable::

    >>> import forgi.utilities.stuff as fus
    >>> fus.dotbracket_to_pairtable('((..))..((..))')
    [14, 6, 5, 0, 0, 2, 1, 0, 0, 14, 13, 0, 0, 10, 9]

The pair-table list has the number of nucleotides as the first element.
Every subsequent entry at position `i` contains the pairing partner of
the `i`'th nucleotide. If the value at position `i` is `0`, then nucleotide
`i` is unpaired.

Pair-table to dot-bracket
-------------------------

Converting a pair-table to a dot-bracket representation::

    >>> import forgi.utilities.stuff as fus
    >>> fus.pairtable_to_dotbracket([14, 6, 5, 0, 0, 2, 1, 0, 0, 14, 13, 0, 0, 10, 9])
    '((..))..((..))'
