#!/usr/bin/python

from __future__ import print_function
from __future__ import division

import sys
import os.path
import subprocess as sp
import logging
import warnings
import argparse
from collections import defaultdict
import numpy as np

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.commandline_utils as fuc
import forgi.utilities.debug as fud
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.model.similarity as ftms
import forgi.threedee.visual.pymol as ftvp
from forgi.utilities.stuff import make_temp_directory

log = logging.getLogger(__name__)


def get_parser():
    parser = fuc.get_rna_input_parser('Visualize a 3D structure in pymol.',
                                      '+', '3d')

    parser.add_argument('--thin-cylinders', default=None,
                        help='Make coarse_grain RNA thinner')
    parser.add_argument('--virtual-atoms', default=False, action='store_true',
                        help='Display the virtual atoms')
    parser.add_argument('--virtual-residues', default=False,
                        action='store_true', help='Display the virtual residues as spheres')
    parser.add_argument('--virtual-stems', default="",
                        type=str, help='Display the virtual stems at broken MLS (in this color)')
    parser.add_argument('--three-points', default=False,
                        action='store_true', help='Display the virtual 2 points as spheres')
    parser.add_argument('--only-elements', dest='only_elements', default=None,
                        help='Display only these elements, separated by commas')
    parser.add_argument('--no-loops', action="store_true",
                        help="Don't display the coarse-grain hairpin loops")
    parser.add_argument('--longrange', default=False,
                        action='store_true', help="Display long-range interactions")
    parser.add_argument('--stem-color', default='green',
                        help='The default stem color in coarse-grain drawings')
    parser.add_argument('--multiloop-color', default='red',
                        help='The default multiloop color in coarse-grain drawings')
    parser.add_argument('-x', '--text', default=False,
                        action='store_true',
                        help="Add labels indicating the element names to the figure.")
    parser.add_argument('--labels', type=str,
                        help="Add labels to elements. Expects a comma seperated "
                        "string of element:label, like 'm0:LookHere' to"
                        " display 'LookHere' at the center of 'm0'.")
    parser.add_argument('--encompassing-stems', default=False,
                        action='store_true',
                        help=argparse.SUPPRESS)  # 'Show the big stems that encompass the colinear ones.')
    parser.add_argument('--sidechain-atoms', default=False, action='store_true',
                        help='Include the sidechain atoms. Automatically enables --virtual-atoms')
    parser.add_argument('--basis', default=False, action='store_true',
                        help=argparse.SUPPRESS)  # Display the second basis vector of each element's coordinate system at its origin")
    parser.add_argument('--rainbow', default=False, action='store_true',
                        help='Color virtual atoms (if displayed) depending on the nucleotide position.')
    parser.add_argument('--element-colors', default="",
                        help='Specify a color for certain elements '
                        '(comma-separated element names or element_name:color)'
                        'Example: "m1,m2" (makes m1, m2 purple) or '
                        '"m0:red,m1:green,s0:AA11GG,default:black", '
                        'where "AA11GG" is a hex value and "default" addresses all other elements.'
                        ' Warning: colors not understood by PYMOL will be interpreted as black.'
                        )
    parser.add_argument('--align', default=False, action='store_true',
                        help='Align the RNAs (if multiple are provided)')
    parser.add_argument('-o', '--output', default=None, help="Create a picture of the scene and exit",
                        type=str)
    parser.add_argument('--batch', default=False, action='store_true',
                        help='Start pymol in batch mode')
    parser.add_argument('--pymol-file', type=str,
                        help=argparse.SUPPRESS)  # Store the PYMOL file under this name. WARNING: Do not use .pml as file ending!!!
    parser.add_argument('--only-first-bg', action="store_true", help=argparse.SUPPRESS)

    return parser


def pymol_printer_from_args(args):
    pp = ftvp.PymolPrinter()
    if args.thin_cylinders:
        pp.cylinder_width = 0.5
        pp.show_twists = False
    pp.display_virtual_residues = args.virtual_residues    
    pp.display_3_points = args.three_points
    pp.virtual_atoms = args.virtual_atoms

    if args.only_elements is not None:
        pp.only_elements = args.only_elements.split(',')
    pp.add_loops = not args.no_loops
    pp.add_longrange = args.longrange
    pp.stem_color = args.stem_color
    pp.multiloop_color = args.multiloop_color
    pp.print_text = args.text
    pp.encompassing_stems = args.encompassing_stems
    pp.sidechain_atoms = args.sidechain_atoms
    pp.basis = args.basis
    pp.rainbow = args.rainbow
    pp.plot_virtual_stems = args.virtual_stems
    if args.virtual_stems:
        if args.virtual_stems not in ftvp.NAMED_COLORS:  # A hex value
            try:
                color = tuple(
                    int(args.virtual_stems[i:i + 2], 16) / 255 for i in (0, 2, 4))
            except:
                raise ValueError("Color value '{}' not understood. "
                                 "Either provide a HEX value or "
                                 "one of {}".format(args.virtual_stems, ",".join(ftvp.NAMED_COLORS.keys())))
            pp.plot_virtual_stems = color
    if args.element_colors:
        directives = args.element_colors.split(",")
        elem_colors = {}
        for directive in directives:
            elem, _, color = directive.partition(":")
            log.debug("Element %s: %s", elem, color)
            if not color:
                color = "purple"
            if color not in ftvp.NAMED_COLORS:  # A hex value
                try:
                    color = tuple(
                        int(color[i:i + 2], 16) / 255 for i in (0, 2, 4))
                except:
                    raise ValueError("Color value '{}' not understood. "
                                     "Either provide a HEX value or "
                                     "one of {}".format(color, ",".join(ftvp.NAMED_COLORS.keys())))
            if elem.lower() == "default":
                default_color = color  # Color is changed later
                elem_colors = defaultdict(lambda: default_color, elem_colors)
            else:
                elem_colors[elem] = color
            log.debug("Element %s", elem_colors)

        pp.element_specific_colors = elem_colors
    return pp


def align_rnas(rnas):
    crds0 = rnas[0].get_ordered_virtual_residue_poss()
    centroid0 = ftuv.get_vector_centroid(crds0)
    print(centroid0)
    rnas[0].rotate_translate(centroid0, ftuv.identity_matrix)
    crds0-=centroid0
    assert  ftuv.magnitude(ftuv.get_vector_centroid(crds0))<10**-5, ftuv.magnitude(ftuv.get_vector_centroid(crds0))
    assert  ftuv.magnitude(ftuv.get_vector_centroid(rnas[0].get_ordered_virtual_residue_poss()))<10**-5, ftuv.get_vector_centroid(rnas[0].get_ordered_virtual_residue_poss())
    for rna in rnas[1:]:
        crds1 = rna.get_ordered_virtual_residue_poss()
        centroid1 = ftuv.get_vector_centroid(crds1)
        crds1-=centroid1
        rot_mat = ftms.optimal_superposition(crds0, crds1)
        rna.rotate_translate(centroid1, rot_mat)
        assert  ftuv.magnitude(ftuv.get_vector_centroid(crds1))<10**-5, ftuv.magnitude(ftuv.get_vector_centroid(crds1))
        assert  ftuv.magnitude(ftuv.get_vector_centroid(rna.get_ordered_virtual_residue_poss()))<10**-5, ftuv.magnitude(ftuv.get_vector_centroid(rna.get_ordered_virtual_residue_poss()))


def main(args):
    rnas = fuc.cgs_from_args(args, '+', '3d')
    pp = pymol_printer_from_args(args)

    if args.align:
        print("Aligning RNAs")
        align_rnas(rnas)
    if args.labels:
        label_list = args.labels.split(",")
        labels = {}
        for label in label_list:
            if not label:
                continue
            try:
                elem, lab = label.split(':')
            except ValueError:
                raise ValueError(
                    "Please specify --labels with as list of colon-seperated tuples. Found invalid entry {}.".format(repr(label)))
            labels[elem] = lab
        if not pp.print_text:
            labels = defaultdict(lambda: "", labels)
            pp.print_text = True
    else:
        labels = {}

    color_modifier = 1.0
    log.info("Visualizing {} rnas".format(len(rnas)))
    plot_bg = True
    for rna in rnas:
        pp.add_cg(rna, labels, color_modifier, plot_core_bulge_graph=plot_bg)
        if args.only_first_bg:
            plot_bg=False
        #color_modifier *= 0.7

    with make_temp_directory() as tmpdir:
        # The file describing the cg-structure as cylinders
        if args.pymol_file:
            stru_filename = args.pymol_file
        else:
            stru_filename = os.path.join(tmpdir, "structure")
        with open(stru_filename, "w") as f:
            f.write(pp.pymol_string())

        pdb_fns = []
        selections = ""
        group_selections=""
        for i, rna in enumerate(rnas):
            sel_names=[]
            if rna.chains:
                obj_name = "pdb{}_{}".format(i, rna.name.replace("-", "_"))
                fn = os.path.join(tmpdir, obj_name + ".cif")
                pdb_fns.append(fn)
                ftup.output_multiple_chains(rna.chains.values(), fn, "cif")
                for d in rna.defines:
                    resids = list(
                        rna.define_residue_num_iterator(d, seq_ids=True))
                    if resids:
                        chains = {r.chain for r in resids}
                        sel = []
                        for c in chains:
                            sel.append("( %{} and chain {} and resi {}) ".format(
                                obj_name, c, "+".join(map(str, (r.resid[1] for r in resids)))))
                        sel_name = d + "_" + obj_name
                        selections += "select {}, {}\n".format(sel_name, " or ".join(sel))
                        sel_names.append(sel_name)
                group_selections+=("cmd.group('sel_{}', '{}')\n".format(obj_name, " ".join(sel_names)))

        pymol_cmd = 'hide all\n'
        pymol_cmd += 'show cartoon, all\n'
        pymol_cmd += 'set cartoon_ring_mode\n'
        pymol_cmd += 'set cartoon_tube_radius, .3\n'
        if args.only_elements is not None:
            pymol_cmd += "hide all\n"

            for constraint in args.only_elements.split(','):
                color = pp.get_element_color(constraint)
                for i, rna in enumerate(rnas):
                    for r in rna.define_residue_num_iterator(constraint, seq_ids=True):
                        pymol_cmd += "show sticks, resi {}\n".format(r[1])
                        pymol_cmd += "color {}, resi {}\n".format(color, r[1])

        pymol_cmd += 'run %s\n' % (stru_filename)
        pymol_cmd += 'bg white\n'
        pymol_cmd += 'clip slab, 10000\n'
        #pymol_cmd += 'orient\n'
        pymol_cmd += selections
        pymol_cmd += group_selections

        if args.output is not None:
            pymol_cmd += 'ray\n'
            pymol_cmd += 'png %s\n' % (args.output)
            #pymol_cmd += 'quit\n'
        pml_filename = os.path.join(tmpdir, "command.pml")
        with open(pml_filename, "w") as f1:
            f1.write(pymol_cmd)
        if args.batch:
            p = sp.Popen(['pymol', '-cq'] + pdb_fns +
                         [pml_filename], stdout=sp.PIPE, stderr=sp.PIPE)
        else:
            p = sp.Popen(['pymol'] + pdb_fns + [pml_filename],
                         stdout=sp.PIPE, stderr=sp.PIPE)
        log.info("Now opening pymol")
        out, err = p.communicate()
        log.info("Out=\n%s", out)
        log.info("Errt=\n%s", err)

parser = get_parser()
if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
