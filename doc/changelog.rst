Changes in upcoming version
===========================
*  `average_stem_vres_positions` are back with recalculated values
*  In `forgi.graph.bulge_graph` the function `BulgeGraph.get_connected_nucleotides` no longer sorts the 
   output nucleotides. Now this function depends on the order of stem1 and stem2 and can thus be used 
   to determine the direction of a bulge. This is used in the new member function `get_link_direction`.


Changes in Version 0.3
======================

*  CoarseGrainRNA now has a member `cg.virtual_atoms()` which is used for caching of virtual atom positions.
   `forgi.threedee.utilities.graph_pdb.virtual_atoms()` now only calculates atom positions when they are needed.
   The two changes together lead to a great speed improvement.
*  The subpackage `forgi.visual` was started for easy visualization of RNA using fornac or matplotlib.
   This subpackage is in an early development stage and will be improved in future versions.
*  The subpackage forgi.projection was added to work with projections of CoarseGrainedRNA objects onto a 2D plane.
*  Now `forgi.threedee.model.average_atom_positions` is used for all virtual atom calculations 
   while `average_stem_vres_positions` has been removed. This leads to more consistent virtual atom calculations.
   Further more the values in `average_atom_positions` have been recalculated.
*  In `forgi.threedee.model.rmsd`, the functions `centered_rmsd` and `centered_drmsd` have been 
   deprecated and deleted respectively. Use `rmsd(cg1,cg2)` for a centered RMSD. This removes code duplication.
*  In `forgi.threedee.model.comparison` a ConfusionMatrix class was introduced for speedup with 
   repeated comparisons to the same reference.
*  Several smaller changes and improvements
