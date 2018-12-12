from __future__ import print_function
from __future__ import division
from builtins import zip
from past.builtins import basestring

import Bio
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.pdb as ftup
import math
import matplotlib.pyplot as plt
import numpy as np
import logging

log = logging.getLogger(__name__)

def circles(x, y, s, c='b', ax=None, vmin=None, vmax=None,labels=[], **kwargs):
    """
    Make a scatter of circles plot of x vs y, where x and y are sequence
    like objects of the same lengths. The size of circles are in data scale.

    :param x,y: scalar or array_like, shape (n, )
                Input data
    :param s: scalar or array_like, shape (n, )
              Radius of circle in data scale (ie. in data unit)
    :param c: color or sequence of color, optional, default : 'b'
              `c` can be a single color format string, or a sequence of color
              specifications of length `N`, or a sequence of `N` numbers to be
              mapped to colors using the `cmap` and `norm` specified via kwargs.
              Note that `c` should not be a single numeric RGB or
              RGBA sequence because that is indistinguishable from an array of
              values to be colormapped.  `c` can be a 2-D array in which the
              rows are RGB or RGBA, however.
    :param ax: Axes object, optional, default: None
               Parent axes of the plot. It uses gca() if not specified.
    :param vmin, vmax: scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.  (Note if you pass a `norm` instance, your
        settings for `vmin` and `vmax` will be ignored.)
    :param kwargs: `~matplotlib.collections.Collection` properties
        eg. alpha, edgecolors, facecolors, linewidths, linestyles, norm, cmap

    :returns: paths : `~matplotlib.collections.PathCollection`


    Examples ::

        a = np.arange(11)
        circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')

    This code by ??? is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)

    """
    from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection
    #import matplotlib.colors as colors

    if ax is None:
        ax = plt.gca()

    if isinstance(c, basestring):
        color = c     # ie. use colors.colorConverter.to_rgba_array(c)
    else:
        color = None  # use cmap, norm after collection is created
    kwargs.update(color=color)

    if np.isscalar(x):
        patches = [Circle((x, y), s), ]
    elif np.isscalar(s):
        patches = [Circle((x_, y_), s) for x_, y_ in zip(x, y)]
    else:
        patches = [Circle((x_, y_), s_) for x_, y_, s_ in zip(x, y, s)]

    collection = PatchCollection(patches, **kwargs)

    if color is None:
        collection.set_array(np.asarray(c))
        if vmin is not None or vmax is not None:
            collection.set_clim(vmin, vmax)

    ax.add_collection(collection)
    ax.autoscale_view()
    return collection


def plot_rna(cg, ax=None, offset=(0, 0), text_kwargs={}, color=True):
    '''
    Plot an RNA structure given a set of nucleotide coordinates

    :param cg: A forgi.threedee.model.coarse_grain.CoarseGrainRNA structure
    :param ax: A matplotlib plotting area
    :param offset: Offset the plot by these coordinates. If a simple True is passed in, then
                   offset by the current width of the plot
    :param text_kwargs: keyword arguments passed to matplotlib.pyplot.annotate.
    :return: (ax, coords) The axes and the coordinates for each nucleotide
    '''
    log.info("Starting to plot RNA...")
    import RNA

    RNA.cvar.rna_plot_type = 1

    coords = []
    #colors = []
    #circles = []

    bp_string = cg.to_dotbracket_string()
    # get the type of element of each nucleotide
    el_string = cg.to_element_string()
    # i.e. eeesssshhhhsssseeee
    el_to_color = {'f': 'orange',
                   't': 'orange',
                   's': 'green',
                   'h': 'blue',
                   'i': 'yellow',
                   'm': 'red'}

    if ax is None:
        ax = plt.gca()

    if offset is None:
        offset = (0, 0)
    elif offset is True:
        offset = (ax.get_xlim()[1], ax.get_ylim()[1])
    else:
        pass

    vrna_coords = RNA.get_xy_coordinates(bp_string)
    for i, _ in enumerate(bp_string):
        coord = (offset[0] + vrna_coords.get(i).X,
                 offset[1] + vrna_coords.get(i).Y)
        coords.append(coord)
        #colors += [el_to_color[el_string[i-1]]]
        if color:
            circle = plt.Circle((coord[0], coord[1]),
                            color=el_to_color[el_string[i]])
        else:
            circle = plt.Circle((coord[0], coord[1]),
                                edgecolor="black", fill=False)

        ax.add_artist(circle)
        if cg.seq:
            if "fontweight" not in text_kwargs:
                text_kwargs["fontweight"]="bold"
            ax.annotate(cg.seq[i+1],xy=coord, ha="center", va="center", **text_kwargs )

    coords = np.array(coords)
    datalim = ((min(list(coords[:, 0]) + [ax.get_xlim()[0]]),
                min(list(coords[:, 1]) + [ax.get_ylim()[0]])),
               (max(list(coords[:, 0]) + [ax.get_xlim()[1]]),
                max(list(coords[:, 1]) + [ax.get_ylim()[1]])))

    '''
    min_coord = min(datalim[0][0], datalim[0][1])
    max_coord = max(datalim[1][0], datalim[1][1])
    datalim = ((min_coord, min_coord), (max_coord, max_coord))

    print "min_coord:", min_coord
    print "max_coord:", max_coord
    print "datalime:", datalim
    '''

    width = datalim[1][0] - datalim[0][0]
    height = datalim[1][1] - datalim[0][1]

    #ax.set_aspect(width / height)
    ax.set_aspect('equal', 'datalim')
    ax.update_datalim(datalim)
    ax.autoscale_view()

    return (ax, coords)


def plot_pdb(filename, ax=None):
    """
    Plot a pdb file.

    :param structure: A Bio.PDB.Structure
    :return: An Axes object (ax)
    """

    structure = Bio.PDB.PDBParser().get_structure('blah', filename)
    model = list(structure)[0]
    ax = None
    chain_coords = {}
    cgs = {}

    import collections as col

    # store a list of RNA nucleotides that each protein interacts with
    protein_interactions = col.defaultdict(set)
    protein_circles = []

    for chain in model:
        # iterate over RNAs
        if ftup.is_rna(chain):
            # convert to cg and store so that we can convert pdb nucleotide ids to
            # secondary structure indexes later
            cg = ftmc.from_pdb(filename, chain_id=chain.id)
            cgs[chain.id] = cg

            # plot the structure and store the coordinates
            (ax, coords) = plot_rna(cg, offset=True, ax=ax)
            chain_coords[chain.id] = coords

    for (a1, a2) in ftup.interchain_contacts(structure):
        # iterate over all the interactions in order to find out which
        # nucleotides this protein interacts with
        chain1 = a1.parent.parent
        chain2 = a2.parent.parent

        if ftup.is_protein(chain1) and ftup.is_rna(chain2):
            # collect all the RNA nucleotides that a protein interacts with
            sid = cgs[chain2.id].seq_ids.index(a2.parent.id)
            protein_interactions[chain1.id].add((chain2.id, sid))

        if ftup.is_rna(chain1) and ftup.is_rna(chain2):
            sid1 = cgs[chain1.id].seq_ids.index(a1.parent.id)
            sid2 = cgs[chain2.id].seq_ids.index(a2.parent.id)

            coord1 = chain_coords[chain1.id][sid1]
            coord2 = chain_coords[chain2.id][sid2]

            ax.plot([coord1[0], coord2[0]], [coord1[1], coord2[1]],
                    'k-', alpha=0.5)

    for chain in model:
        # draw each protein and the links that it has to other nucleotides
        if ftup.is_protein(chain):
            # the protein will be positioned at the centroid of the nucleotides
            # that it interacts with
            interacting_coords = [np.array(chain_coords[chain_id][nuc_num])
                                  for (chain_id, nuc_num) in protein_interactions[chain.id]]

            centroid = np.sum(interacting_coords, axis=0) / \
                len(interacting_coords)

            # the size of the circle representing it will be proportional to its
            # length (in nucleotides)
            radius = 2 * math.sqrt(len(chain.get_list()))
            protein_circles += [[centroid[0], centroid[1], radius]]

            # draw all of the interactions as lines
            for coord in interacting_coords:
                ax.plot([coord[0], centroid[0]], [
                        coord[1], centroid[1]], 'k-', alpha=0.5)

    protein_circles = np.array(protein_circles)
    if len(protein_circles) > 0:
        circles(protein_circles[:, 0], protein_circles[:, 1],
                protein_circles[:, 2], 'grey', alpha=0.5)

    # plt.axis('off')

    pass
