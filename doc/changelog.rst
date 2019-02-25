Changelog
=========

Changes in Version 2.0
----------------------

Version 2.0 brings a cleaned-up API to forgi. The BulgeGraph object should now be treated as immuteable on the level of defines and edges. Furthermore, the **new Sequence class** deals with many quirks of the PDB format, including tuple-indices and missing and modified residues.
We now **fully support the MMCIF fileformat**.


Some of the biggest changes:

*  All RNA file-formats can be **loaded with the single function `forgi.load_rna`**, which returns a list of connected components.
*  The old functions for creating the RNA from a specific file type are now classmethods of the BulgeGraph class. `BulgeGraph.__init__` should no longer be used.
* The API of the Bulgegraph object was cleaned up, including the following changes:

  *  Functions modifying the defines and edges of the BulgeGraph (like `add_node`) have been moved to the private module forgi.graph._graph_construction, as they are only used during creation of the BulgeGraph object, including `dissolve_length_one_stems`. Instead, `length_one_stem_basepairs` was added.
  *  To get a modified copy of a BulgeGraph, use the member methods of the GraphTransformer, accessed with `.transformed`. These methods are implemented in forgi/graph/transform_graphs.py and more methods will be added in the future.
  *  Some functions that are probably not useful for client code have been made private or removed, including:

    *  `all_connections` -> `_all_connections`
    *  `compare_stems`, `compare_*`: removed
    *  `find_external_loops`, use `get_domains`
    *  `flanking_nucleotides`, use `define_a`.
    *  `get_any_sides`, `get_sides_plus` removed
    *  Functions only called by `to_bg_string` are now private 
       (`get_*_string`, e.g. `_get_connect_str`).
    *  `get_stem_direction` removed
    *  `get_vertex` removed
    *  `nd_define_iterator` removed
    *  `subseq_with_cutpoints` removed, use Sequence class's `__getitem__`
    
  *  `get_node_from_residue_num` is deprecated and shuld be replaced by `get_elem`
  *  Better support for `forgi.graph.residue.RESID` as indices.
  *  New properties/ methods `junctions`, `rods`, `pseudoknotted_basepairs`, `seq_length`
  
*  A **new sequence class**, which allows for indexing with integers (1-based, like the defines) and `forgi.graph.residue.RESID` tuples. Conversion between these two indices is now also delegated to this class. `__getitem__` is implemented to properly allow for slices with cutpoints.
*  Some changes were made to the CoarseGrainRNA class:
  
  *  **Virtual residues are now stored for loops** (at the C1' atom), if the RNA was loaded from a PDB or MMCIF file.
  *  Experimental functions are now made private.
  *  New function `get_incomplete_elements` for elements with missing residues.
  *  New function `rotate_translate`
  
*  PDB files can now be read, even if neither DSSR nor MC-Annotate are installed. However, the accurracy of the secondary structure may suffer in this case.
* `forgi.threedee.model.similarity.cg_rmsd` now takes virtual residues of loops into account.
* The module `forgi.graph.numbered_dotbracket` allows you to modify the dotbracket structure while keeping track of the mapping between brackets and residues.
 
Changes in Version 1.1
----------------------

*  We now support **DSSR** as an alternative to **MC-Annotate**.
   See the :ref:`installation` for more information
*  RESID for pdb-derived residue ids in now moved into a seperate module.
*  `describe_cg.py` now adds the number of the multiloop in the sequence to the generated csv.
*  **Bugfixes**:

   *  python3 compatibility of `visualize_cg.py`
   *  windows-compatibility of `visualize_cg.py`
   *  fix count of multiloop-types in `describe_cg.py`
   *  persist cofold-breakpoints in cg-files, even if no sequence was given.
   *  fix bugs in parsing of MC-Annotate output. Forgi should now be able to load most pdb-files.
   *  properly propagate the `--keep-length-one-stems` option.

Changes in Version 1.0
----------------------

With version 1.0, we introduce support for **cofold** and multifold structures.
Further more, we are finally compatible with **Python 3** and **Python 2**


In the following list of the biggest changes, the most important
changes are highlighted in bold.

*  **Python 3 compatibility:** forgi is now compatible with python 2 and 3.
   In paticular, we have automatic tests for python 3.5, 3.6 and 2.7
*  Logging support. We now use the `logging` module instead of print statements.
*  In `forgi.graph.bulge_graph`:

    *  **Cofold Structures** A BulgeGraph can now hold multiple RNA chains,
        if they are all connected by basepairs.
    *  bg.seq_ids now holds named tuples (chain, resid)
    *  New Errors (`GraphConstructionError` and `GraphIntegrityError`) as
       specializations of ValueError are introduced and sometimes raised
       by BulgeGraph instances.
    *  Fasta sequences containing 'T' characters are now allowed (and converted to 'U')
    *  **Sequence class:** Introduced `Sequence` class, a string subclass with 1-based indexing and support for '&' (cofold structures)
    *  `define_range_iterator` no longer accepts the `seq_ids` keyword.
    *  The Bulge-graph now has a new method `define_a`. It returns the define with
       adjacent nucleotides (if present) and also works for 0-nucleotide multiloops.
    *  **Consistent element numbering:**

        *  During graph creation, multiloops are now numbered from m0, m1, ... according
           to their position along the backbone.
        *  The numbering of elements during bulge-graph construction is now
           consistent independent of the order of dictionaries
           (which was randomized in python 3 and is ordered in python 3.6)
    *  **Looking at multiloops as a whole:**
        Multiloops continue to consist of independent multiloop segments.
        However, we additionally introduced methods to look at multiloops as a whole.

        *  `find_mlonly_multiloops` finds out, which multiloop segment belong
            together in a junction
        *  `describe_multiloop`: reports whether the junctions found with the
            previous method belong to a pseudoknot, a multiloop or an exterior loop.
    *  `get_angle_type` now supports the `allow_broken` keyword to return an
       angle type instread of None in case of Multiloop segments that are not
       part of the Minimum Spanning Tree
    *  `bg.iter_elements_along_backbone` introduced.
*  `forgi.threedee.utilities.average_atom_positions` was converted to a JSON
    data file, which is read by `forgi.threedee.utilities.graph_pdb`. This saves
    time upon import. In the future, average_atom_positions might be entirely removed.
*  New module `forgi.threedee.model.linecloud` to store the coordinate and twist vectors.
   This module's classes implement a Mapping interface but internally use a numpy array,
   which allows for speedup of many operations.
*  In `forgi.graph.bulge_graph`:

    *  Exceptions `CgConstructionErrorCgConstructionError` and `CgIntegrityError` as
       subclasses of the new bulge graph exceptions
    *  **Load all chains from a PDB file**: The module level function
       `connected_cgs_from_pdb` loads all RNA chains from a PDB file
       and greates a cg object for each connected component.
    *  Bugfix: Supplying the secondary structure during loading of PDBs
       now works again.
    *  Convenience function `cg.get_stats` as a wrapper around `get_bulge_abgle_stats`,
       `get_loop_stats` and `get_stem_stats` respectively.
    *  Code for steric value and stacking detection is EXPERIMENTAL!
*  Modules `ftm.ensemble`, `ftm.ensemble2` and `ftu.dssr` are EXPERIMENTAL
*  Faster RMSD by using QCOrot algorithm with Cython instead of Kabsch algorithm
*  AngleStat class now supports unary minus operator.
*  Recalculated `average_stem_vres_positions` from the NR-list.
*  In `forgi.threedee.utilities.graph_pdb`: virtual stats and sum of angle stats
    to describe orientation of non-adjacent stems in multiloops.
*  **Modified Residues**: Better treatement of modified residues.
   We try to query PDBeChem to replace the modified residue with the unmodified parent.




Changes in Version 0.4
----------------------

*  In `forgi.graph.bulge_graph`:

   *  Speed improvement: Basepair distances between elements are cached.
   *  The Bulge-graph object and file format supports arbitrary key-value pairs in the `info` dictionary.
   *  `BulgeGraph.get_connected_nucleotides` no longer sorts the
      output nucleotides. Now this function depends on the order of stem1 and stem2 and can thus be used
      to determine the direction of a bulge. This is used in the new member function `get_link_direction`.
   *  Added function `BulgeGraph.get_domains`, to return a lists of pseudoknots, multiloops and rods.
      The interface of this function might change in the future.
   *  Merged pull-request by tcarlile for `forgi.graph.bulge_graph`:

      *  `BulgeGraph.get_stem_edge(self, stem, pos)`: Returns 0 if pos on 5' edge of stem, returns 1 if pos on 3' edge of stem.
      *  `BulgeGraph.shortest_path(self, e1, e2)`:    Returns a list of the nodes along the shortest path between e1, and e2.

*  Restructured forgi.threedee.model.comparison and forgi.threedee.utilities.rmsd into
   `forgi.threedee.model.similarity` and `forgi.threedee.model.descriptors`
   The `similarity` module contains all functions for the comparison of two point
   clouds or two cg structures.
   The `descriptors` module contains functions for describing a single point cloud,
   such as the radius of gyration or new functions for the gyration tensor.
*  `average_stem_vres_positions` are back with recalculated values
*  Changes in `forgi.threedee.model.coarse_grain` to the `CoarseGrainRNA` object:

   *  In the `self.coords` dictionary, the start and end coordinate are now in a consistent order.
   *  Call new member function `self.add_all_virtual_residues` instead of `forgi.threedee.utilities.graph_pdb.add_virtual_residues`
   *  Coordinates of interior loops and multiloop segmentsa are no longer stored in the cg-files,
      as they can be deduced from stem coordinates.

      * New member function `self.add_bulge_coords_from_stems` is provided instead
        of function `forgi.threedee.utilities.graph_pdb.add_bulge_information_from_pdb_chain`

   *  Member function `self.get_virtual_residue(pos)` is provided for easier access than directly via `self.v3dposs`.
      For single stranded regions, virtual residue positions along the direct line of the coarse
      grain element can be estimated optionally.
      Virtual residues are cached and the cache is automatically cleared when
      the coordinates or twists of the coarse grained RNA are changed.
   * Functions for creating coordinate arrays for the structure

      * `self.get_ordered_stem_poss` for the start and end coordinates of stems.
      * `self.get_ordered_virtual_residue_poss` for all virtual residue coordinates of stems.
        Replaces `forgi.threedee.utilities.graph_pdb.bg_virtual_residues`
      * `self.get_poss_for_domain` for coordinates only for certain coarse grained elements.

   *  Removed the addition of a pseudo-vector to loop stats in `get_loop_stats`, which was used to avoid zero-length bulges.
      Now 0-length bulges are possible. This makes saving and loading stats consistent.
   *  `self.get_coordinates_array` now returns a 2D `nx3 numpy array` holding n coordinate entries.
      You can use numpy's `.flatten()` to generate a 1D array. If you want to load a 1D flat coordinate array `a`, use
      `self.load_coordinates_array(a.reshape((-1,3))`
   *  Overrides the newly introduced method `sorted_edges_for_mst` from BulgeGraph.
      Now elements that have no `sampled` entry are broken preferedly.
      This should ensure that the minimal spanning tree is the same after saving and loading an
      RNA generated with the program Ernwin to/from a file.
   *  `self.coords_to_directions` and `coords_from_directions`:
      Export the coordinates as an array of directions (end-start).
      This array has only 1 entry per coarse grained element.

*  In `forgi.threedee.model.stats`: Added class for clustered angle stats.
*  Changes in `forgi.projection.hausdorff`.
*  Changes in `forgi.projection.projection2d`

   *  Faster rotation and rasterization.
   *  selected virtual residues can be included in the projection

*  In `forgi.threedee.utilities.graph_pdb`: Added functions `get_basepair_center` and `get_basepair_plane`.
   This will be used in the future for stacking detection.

Changes in Version 0.3
----------------------

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
