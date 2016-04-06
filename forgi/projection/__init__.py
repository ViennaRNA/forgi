"""
The package forgi.projection is a subpackage of forgi responsible for creating 2D projections of 
an RNA 3D structure. 

The class .projection2d.Projection2D is initialized from a 
forgi.threedee.model.coarse_grain.CoarseGrainRNA object. 
It contains methods for plotting and analyzing the projection.

The module forgi.projection.hausdorff consists of two parts.
On the one hand there are very general functions for calculating the Hausdorff distance
between two images represented as boolean matrices. On the other hand there are heuristic 
functions built on top of the before mentioned functions, responsible for calculating the optimal 
alignment of a CoarseGrainRNA and an image, so that the Hausdorff distance between the projection
of the RNA and the image is minimal.
"""
