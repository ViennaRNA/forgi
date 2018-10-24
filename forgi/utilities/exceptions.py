"""
A module where all costum exceptions should be defines.

We put exceptions into this seperate module, so they can be
imported everywhere without risking circular imports
"""


class GraphConstructionError(ValueError):
    """
    Exceptions raised if the BulgeGraph could be constructed from the given input.
    """
    pass


class GraphIntegrityError(ValueError):
    """
    Exception raised if a BulgeGraph was found to be in an inconsistent or faulty state.
    """
    pass


class CgConstructionError(GraphConstructionError):
    """
    Exceptions raised if the CoarseGrainRNA could be constructed from the given input.

    Raised for Errors related to the 3D information.
    """
    pass


class CgIntegrityError(GraphIntegrityError):
    """
    Exception raised if a BulgeGraphCoarseGrainRNA was found to be in an inconsistent or faulty state,
    related to the 3D Information
    """
    pass
