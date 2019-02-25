Forgi API documentation
=======================

This is the documentation of the forgi package,
generated from the docstrings.

.. note::
    Some private submodules contain docstrings in a different format than the rest of the code base.
    They have been left out of this documentation. Please see the source code for their documentation.

To get to know forgi, you might want to start with reading the documentation of `forgi.graph.bulge_graph` and `forgi.threedee.models.coarse_grain`.

To create an RNA object, the preferred method is the function `forgi.load_rna` instead of using the constructors in the RNA classes directly.
With standard parameters, this function returns a list of `forgi.threedee.model.coarse_grain.CoarsegrainRNA` objects.

Subpackages
-----------

.. toctree::

    forgi.graph
    forgi.projection
    forgi.threedee
    forgi.utilities
    forgi.visual

Submodules
----------

.. toctree::

   forgi.config

Module contents
---------------

.. automodule:: forgi
    :members:
    :undoc-members:
    :show-inheritance:
