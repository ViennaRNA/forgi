#!/usr/bin/python

from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from builtins import object

import sys

import itertools as it
import math as m
import numpy as np
import uuid
import collections as col
import warnings
import logging

import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.model.stats as ftmstat
#import forgi.threedee.utilities.average_stem_vres_atom_positions as cua
import forgi.utilities.debug as fud
import forgi.threedee.utilities.vector as cuv
import forgi.threedee.utilities.vector as ftuv

import Bio.PDB.Model as bpm
import Bio.PDB.Structure as bps
import Bio.PDB as bp

log = logging.getLogger(__name__)

NAMED_COLORS = {
    'green':       [0.0, 1.0, 0.0],
    'forest':      [0.2, 0.6, 0.2],
    'blue':        [0.0, 0.0, 1.0],
    'red':         [1.0, 0.0, 0.0],
    'light red':   [0.9, 0.6, 0.6],
    'orange':      [1., 165 / 255., 0.],
    'yellow':      [1.0, 1.0, 0.0],
    'purple':      [1.0, 0.0, 1.0],
    'white':       [1.0, 1.0, 1.0],
    'cyan':        [0.0, 1.0, 1.0],
    'magenta':     [249 / 255., 132 / 255., 229 / 255.],
    'light gray':  [.8, .8, .8],
    'dark gray':   [.1, .1, .1],
    'middle gray': [.6, .6, .6],
    'gray':        [.6, .6, .6],
    'black':       [0., 0., 0.]
}


def pymol_color(color, modifier):
    def get_color_vec(color):
        if color in NAMED_COLORS:
            return NAMED_COLORS[color]
        else:
            return [0.0, 0.0, 0.0]
    if not isinstance(color, (list, tuple)):
        color = [str(c * modifier) for c in get_color_vec(color)]
    else:
        color = [str(c) for c in color[:3]]
    return color


class PyMolRNA(object):
    def __init__(self, name, color_modifier=1.0):
        self.name = name.replace("-", "_").replace(".", "_")
        self.segments = []
        self.boxes = []
        self.labels = []
        self.spheres = []
        self.cones = []
        self.color_modifier = color_modifier

    def add_sphere(self, p, color='green', width=0.2, text=""):
        self.spheres += [(np.array(p), color, width, text)]

    def add_segment(self, p, n, color='green', width=0.2, text="", key=''):
        """
        :param p: start coordinates
        :param n: end-coordinates
        :param color: A color. Either a 3-element list with values
                     from 0 to 1 for R, G and B,
                     or a color name (as defined in this modules variable NAMED_COLORS)
        """
        self.segments += [(np.array(p), np.array(n), color, width, text, key)]

    def add_dashed(self, point1, point2, width=0.1, color="purple"):
        '''
        Add a dashed line from point1 to point2.
        '''
        dash_length = width / 2
        gap_length = dash_length * 2
        direction = ftuv.normalize(point2 - point1)

        num_dashes = ftuv.magnitude(
            point2 - point1) / (dash_length + gap_length)
        key = None

        for i in range(int(num_dashes)):
            self.add_segment(point1 + i * (dash_length + gap_length) * direction,
                             point1 + (i * (dash_length + gap_length) +
                                       dash_length) * direction, color,
                             width, "", key=key)

    def add_cone(self, p, n, color='white', width=2.4, text=''):

        cone_extension = 2.
        cyl_vec = cuv.normalize(n - p)
        cyl_len = cuv.magnitude(n - p)

        new_width = width * (cyl_len + cone_extension) / cyl_len

        self.cones += [(np.array(p) - cone_extension * cyl_vec,
                        np.array(n), color, width, text)]
        self.cones += [(np.array(n) + cone_extension * cyl_vec,
                        np.array(p), color, width, text)]

    def pymol_segments_string(self):
        color = 'green'
        width = 0.2
        s = ''

        for seg in self.segments:
            (p, n, color, width, text, key) = seg
            color = pymol_color(color, self.color_modifier)
            s += "\n# {}\n".format(key)
            s += " CYLINDER, %f, %f, %f, %f, %f, %f, " % (p[0], p[1], p[2],
                                                          n[0], n[1], n[2])
            s += "%f, %s, %s," % (width, ", ".join(color),
                                  ", ".join(color)) + '\n'
        return s

    def pymol_spheres_string(self):
        s = ''

        for (p, color, width, text) in self.spheres:
            color = pymol_color(color, self.color_modifier)
            if np.ndim(p) > 0:
                s += "COLOR, %s," % (",  ".join([str(c) for c in color]))
                s += '\n'
                s += "SPHERE, %s, %f," % (", ".join([str(pi) for pi in p]),
                                          width)
                s += '\n'
            else:
                warnings.warn(
                    "p is not iterable! It is {} (for '{}'). IGNORING. ".format(p, text))
        return s

    def pymol_text_string(self, rna_number):
        counter = 0
        s = ''
        pa_s = 'cmd.set("label_size", 20)\n'
        uids = []
        names = []
        for (p, n, color, width, text, key) in self.segments:
            if len(text) == 0:
                continue

            # generate a unique identifier for every object so that other
            # scripts can add others that don't clash
            uid = str(uuid.uuid4()).replace('-', 'x')
            uids += [uid]

            if np.all(n == p):
                pos = n
                axes = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
            else:
                comp1 = cuv.normalize(n - p)
                ncl = cuv.get_non_colinear_unit_vector(comp1)
                comp2 = cuv.normalize(np.cross(ncl, comp1))
                comp3 = cuv.normalize(np.cross(ncl, comp2))
                pos = (p + n) / 2.0 + 3 * comp2
                axes = [list(comp1 * 2), list(comp2 * 2), list(comp3 * 2)]

            #text = "%s: %.1f" % (text, cuv.magnitude(n - p))

            name = "label{}_{}_{}".format(len(uids), rna_number, self.name)
            names.append(name)
            pa_s += 'pa_{} = cmd.pseudoatom("{}", pos={},'.format(uid,
                                                                  name, str(list(pos)))
            pa_s += 'b=1.0, label="{}")\n'.format(text)
            counter += 1

        all_patoms = ["pa_{}".format(ui) for ui in uids]
        pa_s += "cmd.group('{}', '{}')\n".format(
            "labels_{}_{}".format(rna_number, self.name), " ".join(names))
        return pa_s

    def pymol_box_string(self):
        '''
        Pring out the CGO text to describe the boxes.
        '''
        out_str = ''
        for (box, color) in self.boxes:
            uid = str(uuid.uuid4()).replace('-', 'x')
            color = pymol_color(color, self.color_modifier)
            out_str += 'obj%s = [\n' % (uid)
            out_str += "LINEWIDTH, .8, \n"
            out_str += "BEGIN, TRIANGLES, \n"
            out_str += "COLOR, %s," % (",  ".join([str(c) for c in color]))
            out_str += '\n'
            for corner in box:
                out_str += "VERTEX, %f, %f, %f, \n" % (corner[0],
                                                       corner[1],
                                                       corner[2])
            out_str += 'END \n'
            out_str += '] \n'
            out_str += "cmd.load_cgo(obj%s, 'ss%s')\n" % (uid, uid)

        return out_str


class PymolPrinter(object):
    def __init__(self):
        self.display_virtual_residues = False
        self.plot_virtual_stems = False
        self.display_3_points = False
        self.rainbow = False
        self.basis = None
        self.visualize_three_and_five_prime = True
        self.encompassing_stems = False
        self.state = 2
        self.virtual_atoms = False
        self.sidechain_atoms = False
        self.override_color = None
        self.element_specific_colors = None
        self.print_text = True
        self.add_longrange = False
        self.add_loops = True
        self.max_stem_distances = 0
        self.add_letters = False
        self.draw_axes = False
        self.draw_segments = True
        self.movie = False
        self.stem_color = 'green'
        self.multiloop_color = 'red'
        self.prev_obj_name = ''     # The name of the previously created
        # object which needs to be hidden
        # when creating a movie
        self.only_elements = None
        self.cylinder_width = 1.0
        self.show_twists = True
        self.plotters = []
        self.show_bounding_boxes = False
    def add_cg(self, cg, labels, color_modifier=1.0, plot_core_bulge_graph=True):
        """
        :param labels: A dictionary with element names as keys
                       and labels as values.
        """
        rna_plotter = PyMolRNA(cg.name, color_modifier)
        for key in cg.coords.keys():
            if self.only_elements is not None:
                if key not in self.only_elements:
                    continue

            (p, n) = cg.coords[key]
            color = self.get_element_color(key)

            if key[0] == 's':
                try:
                    text = labels[key]
                except KeyError:
                    text = key
                if plot_core_bulge_graph:
                    self.add_stem_like(rna_plotter, cg, text, key, color=color)
                if self.show_bounding_boxes:
                    self.draw_bounding_boxes(rna_plotter, cg, key)
            else:
                if key[0] == 'h':
                    if self.add_loops:
                        try:
                            text = labels[key]
                        except KeyError:
                            text = key + " " + str(cg.get_length(key))
                        if plot_core_bulge_graph:
                            rna_plotter.add_segment(p, n, color, self.cylinder_width,
                                                    text,
                                                    key=key)
                elif key[0] == 'm':
                    twists = cg.get_twists(key)
                    try:
                        text = labels[key]
                    except KeyError:
                        # check if the multiloop is longer than one. If it's not, then
                        # it has an empty define and its length will be 0
                        if len(cg.defines[key]) == 0:
                            text = key + " 0"
                        else:
                            text = key + " " + \
                                str(cg.defines[key][1] -
                                    cg.defines[key][0] + 1)
                    if plot_core_bulge_graph:
                        rna_plotter.add_segment(p, n, color, self.cylinder_width,
                                            text, key=key)
                    vstat_line = cg.infos["vstat_{}".format(key)]
                    if self.plot_virtual_stems and key not in cg.get_mst() and vstat_line:
                        log.info("Plotting virtual stem for %s", key)
                        virtual_stat = ftmstat.AngleStat()
                        vstat_line = vstat_line[0]
                        print(vstat_line)
                        virtual_stat.parse_line(vstat_line)
                        stems= cg.edges[key]
                        fixed_stem_name, _ = sorted(stems, key=cg.buildorder_of)
                        vbulge_vec, vstem_coords0, vstem_vec = ftug._get_vstem_coords(cg, key, fixed_stem_name, virtual_stat)
                        rna_plotter.add_segment(vstem_coords0, vstem_coords0+vstem_vec, self.plot_virtual_stems,
                                                2*self.cylinder_width, "", key="virtual stem {}".format(key))
                        fixed_side = cg._get_sides_plus(fixed_stem_name, key)[1]
                        rna_plotter.add_segment(cg.coords[key][fixed_side], vstem_coords0, self.plot_virtual_stems, 0.7*self.cylinder_width,
                                                "", key="virtual stem {}".format(key))
                    else:
                        log.info("NOT Plotting virtual stem for %s: %s, %s, %s", key, self.plot_virtual_stems, key not in cg.get_mst(), vstat_line )

                elif key[0] in 'ft':
                    try:
                        text = labels[key]
                    except KeyError:
                        text = key + " " + \
                            str(cg.defines[key][1] - cg.defines[key][0] + 1)

                    if self.visualize_three_and_five_prime:
                        if plot_core_bulge_graph:
                            rna_plotter.add_segment(p, n, color, self.cylinder_width,
                                                    text, key=key)
                elif key[0] == "i":
                    try:
                        text = labels[key]
                    except KeyError:
                        text = key
                    if plot_core_bulge_graph:
                        rna_plotter.add_segment(
                                p, n, color, self.cylinder_width, text, key=key)

        if self.display_virtual_residues:
            for i in range(1, cg.seq_length + 1):
                elem =  cg.get_node_from_residue_num(i)
                if not self.only_elements or elem in self.only_elements:
                    try:
                        c = self.element_specific_colors[elem]
                    except (KeyError, TypeError):
                        if elem[0] == "s":
                            c = "cyan"
                        else:
                            c = "magenta"

                    pos = cg.get_virtual_residue(i, True)
                    if i>1 and (i-1) not in cg.backbone_breaks_after:
                        prev_elem =  cg.get_node_from_residue_num(i-1)
                        if elem not in cg.get_mst() or prev_elem not in cg.get_mst():
                            w=0.05
                        else:
                            w=0.2
                        rna_plotter.add_segment(cg.get_virtual_residue(i-1, True), pos, c, w, key="vres{}-{} {}".format(i-1, i, elem))
                    j=cg.pairing_partner(i)
                    if j:
                        rna_plotter.add_segment(cg.get_virtual_residue(j, True), pos, "gray", 0.3, key="vres BP {}-{}-{}".format(j, i, elem))
                    rna_plotter.add_sphere(pos, c, 1.)

        if self.display_3_points:
            for i in range(1, cg.seq_length + 1):
                elem =  cg.get_node_from_residue_num(i)
                if not self.only_elements or elem in self.only_elements:
                    for pos in cg.iter_three_points(i):
                        if cg.get_node_from_residue_num(i)[0] == "s":
                            c = "cyan"
                        else:
                            c = "magenta"
                        rna_plotter.add_sphere(pos, c, 0.4)

        if self.add_longrange:
            for key1 in cg.longrange.keys():
                for key2 in cg.longrange[key1]:
                    if self.only_elements is not None:
                        if key1 not in self.only_elements or key2 not in self.only_elements:
                            continue
                    try:

                        p = cuv.line_segment_distance(cg.coords[key1][0],
                                                      cg.coords[key1][1],
                                                      cg.coords[key2][0],
                                                      cg.coords[key2][1])

                        rna_plotter.add_dashed(p[0], p[1])
                    except:
                        continue

        if self.encompassing_stems:
            self.add_encompassing_cylinders(rna_plotter, cg, 7.)

        if self.max_stem_distances > 0:
            for (s1, s2) in it.permutations(cg.stem_iterator(), r=2):
                (i1, i2) = cuv.line_segment_distance(cg.coords[s1][0],
                                                     cg.coords[s1][1],
                                                     cg.coords[s2][0],
                                                     cg.coords[s2][1])
                if cuv.magnitude(i2 - i1) < self.max_stem_distances:
                    #self.add_segment(i1, i2, 'cyan', 0.3, s1 + " " + s2, key=key)
                    rna_plotter.add_segment(i1, i2, 'cyan', 0.3, key=key)

        if self.virtual_atoms or self.sidechain_atoms:
            cg.add_all_virtual_residues()
            va = ftug.virtual_atoms(cg, sidechain=self.sidechain_atoms)

            atom_width = 0.5
            for i, r in enumerate(sorted(va.keys())):
                for a in va[r].keys():
                    if self.rainbow:
                        import matplotlib
                        matplotlib.use('Agg')
                        import matplotlib.pyplot as plt
                        cmap = plt.get_cmap('gist_rainbow')
                        rna_plotter.add_sphere(va[r][a],
                                               color=cmap(
                                                   i / float(len(va.keys()))),
                                               width=atom_width)
                    else:
                        d = cg.get_node_from_residue_num(r)
                        if d[0] == 's':
                            if a in ftup.nonsidechain_atoms:
                                rna_plotter.add_sphere(
                                    va[r][a], self.stem_color, width=atom_width)
                            else:
                                rna_plotter.add_sphere(
                                    va[r][a], 'forest', width=atom_width)
                        elif d[0] == 'i':
                            rna_plotter.add_sphere(
                                va[r][a], 'yellow', width=atom_width)
                        elif d[0] == 'm':
                            rna_plotter.add_sphere(
                                va[r][a], self.multiloop_color, width=atom_width)
                        elif d[0] == 'h':
                            rna_plotter.add_sphere(
                                va[r][a], 'blue', width=atom_width)

        if self.basis:
            for d in cg.defines.keys():
                origin, basis = ftug.element_coord_system(cg, d)

                rna_plotter.add_segment(
                    origin, origin + 7. * basis[1], 'purple', 0.5, key=key)
        self.plotters.append(rna_plotter)

    def add_stem_like(self, rna_plotter, cg, text, key, color='green', width=2.4):
        if key in cg.twists:
            return self.add_stem_like_core(rna_plotter, cg.coords[key], cg.twists[key],
                                           cg.stem_length(key), text, key, color, width)
        else:
            return self.add_stem_like_core(rna_plotter, cg.coords[key], None,
                                           cg.stem_length(key), text, key, color, width)

    def add_stem_like_core(self, rna_plotter, coords, twists, stem_len, text, key,
                           color='green', width=2.4):
        (p, n) = coords
        width *= self.cylinder_width
        rna_plotter.add_segment(p, n, color, width, text, key=key)

        if self.show_twists:
            rna_plotter.add_cone(p, n, 'white', width, key)
            mult = 8.
            width = .3
            (twist1o, twist2o) = twists

            rna_plotter.add_segment(
                p, p + mult * twist1o, "cyan", width, '', key=key)
            rna_plotter.add_segment(
                n, n + mult * twist2o, "magenta", width, '', key=key)

            for i in range(stem_len):
                res = ftug.virtual_res_3d_pos_core((p, n), twists, i, stem_len)
                (pos, vec_c, vec_l, vec_r) = res
                rna_plotter.add_segment(
                    pos, pos + mult * vec_c, "orange", width, '', key=key)

                if self.add_letters:
                    rna_plotter.labels += [('L', list(pos + mult * vec_l))]
                    rna_plotter.labels += [('R', list(pos + mult * vec_r))]

    def pymol_axis_string(self):
        w = 0.12  # cylinder width
        l = 10.0  # cylinder length
        h = 3.0  # cone hight
        d = w * 2.618  # cone base diameter
        s = ""

        s += "CYLINDER, 0.0, 0.0, 0.0,   %f, 0.0, 0.0, %f," % (l, w)
        s += " 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,"
        s += "CYLINDER, 0.0, 0.0, 0.0, 0.0,   %f, 0.0, %f, " % (l, w)
        s += "0.0, 0.0, 1.0, 0.0, 0.0, 1.0,"
        s += "CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   %f, %f, " % (l, w)
        s += "1.0, 0.0, 0.0, 1.0, 0.0, 0.0,"
        s += "CONE,   %f, 0.0, 0.0, %f, 0.0, 0.0, %f, " % (l, h + l, d)
        s += "0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0,"
        s += "CONE, 0.0, %f, 0.0, 0.0, %f, 0.0, %f, " % (l, h + l, d)
        s += "0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0,"
        s += "CONE, 0.0, 0.0, %f, 0.0, 0.0, %f, %f, " % (l, h + l, d)
        s += "0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0,"

        return s

    def pymol_string(self):
        '''
        Output the contents of this structure into a file that can be passed
        in as a pymol script.
        '''

        s = self.pymol_intro_string()

        for i, plotter in enumerate(self.plotters):
            object_id = "{}_{}".format(i, plotter.name)
            s += "obj{} = [\n".format(object_id)
            if self.draw_segments:
                s += plotter.pymol_segments_string()
            s += plotter.pymol_spheres_string()
            s += ']\n'
            s += self.pymol_load_cgo_string(object_id)

            if self.print_text:
                s += plotter.pymol_text_string(i)
            s += plotter.pymol_box_string()

        if self.draw_axes:
            s += plotter.pymol_axis_string()
        return s

    def pymol_intro_string(self):
        s = "from pymol.cgo import *" + '\n'
        s += "from pymol import cmd" + '\n'
        s += "from pymol.vfont import plain" + '\n'
        return s

    def pymol_load_cgo_string(self, object_id):
        if self.movie:
            s = "cmd.load_cgo(obj%s, 'ss%s', %d)" % (object_id,
                                                     object_id,
                                                     self.state) + '\n'
            self.state += 1
        else:
            s = "cmd.load_cgo(obj%s, 'ss%s')" % (object_id,
                                                 object_id) + '\n'
        return s

    def draw_bounding_boxes(self, rna_plotter, bg, s):
        '''
        Draw bounding boxes for all of the residues encompassed
        by a stem. But only if there is a pdb file handy.

        @param bg: The BulgeGraph
        @param s: The name of the stem
        '''
        if not bg.chains:
            return

        for i in range(bg.stem_length(s)):
            (origin, basis, bb) = ftug.bounding_boxes(bg, s, i)
            for k in range(2):
                (n, x) = bb[k]

                corners = [
                          [n[0], n[1], n[2]],
                          [n[0], n[1], x[2]],

                          [n[0], x[1], n[2]],
                          [n[0], x[1], x[2]],

                          [x[0], n[1], n[2]],
                          [x[0], n[1], x[2]],

                          [x[0], x[1], n[2]],
                          [x[0], x[1], x[2]],

                          [n[0], n[1], n[2]],
                          [x[0], n[1], n[2]],

                          [n[0], x[1], n[2]],
                          [x[0], x[1], n[2]],

                          [n[0], x[1], x[2]],
                          [x[0], x[1], x[2]],

                          [n[0], n[1], x[2]],
                          [x[0], n[1], x[2]],

                          [n[0], n[1], n[2]],
                          [n[0], x[1], n[2]],

                          [x[0], n[1], n[2]],
                          [x[0], x[1], n[2]],

                          [n[0], n[1], x[2]],
                          [n[0], x[1], x[2]],

                          [x[0], n[1], x[2]],
                          [x[0], x[1], x[2]]]

                new_corners = []
                for corner in corners:
                    new_corners += [origin +
                                    cuv.change_basis(np.array(corner),
                                                     cuv.standard_basis,
                                                     basis)]
                corners = np.array(new_corners)

                if k == 0:
                    rna_plotter.boxes += [(corners, 'yellow')]
                    rna_plotter.add_sphere(corners[0], 'yellow', 0.4, '')
                    rna_plotter.add_sphere(corners[7], 'yellow', 0.4, '')
                else:
                    rna_plotter.add_sphere(corners[0], 'purple', 0.4, '')
                    rna_plotter.add_sphere(corners[7], 'purple', 0.4, '')
                    rna_plotter.boxes += [(corners, 'purple')]

    def add_encompassing_cylinders(self, rna_plotter, cg, radius=7.):
        cylinders_to_stems = ftug.get_encompassing_cylinders(cg, radius)

        for stems in cylinders_to_stems.values():
            print("stems:", stems)

            points = []
            for s in stems:
                points += [cg.coords[s][0], cg.coords[s][1]]

            # create the linear regression
            data = np.array(points)
            datamean = data.mean(axis=0)

            uu, dd, vv = np.linalg.svd(data - datamean)

            furthest = max([ftuv.magnitude(d) for d in (data - datamean)])

            start_point = -furthest * vv[0] + datamean
            end_point = furthest * vv[0] + datamean

            rna_plotter.add_segment(
                start_point, end_point, 'white', width=4, text='', key='')

    def get_element_color(self, elem_name):
        '''
        Get the color for this element. The color is determined by the name
        of the element.

        @param elem_name: The name of the element.
        @return: A string with a color name
        '''
        if self.element_specific_colors is not None:
            try:
                return self.element_specific_colors[elem_name]
            except KeyError:
                pass

        if elem_name[0] == 's':
            return self.stem_color
        elif elem_name[0] == 'i':
            return 'yellow'
        elif elem_name[0] == 'm':
            return self.multiloop_color
        elif elem_name[0] == 'h':
            return 'blue'
        elif elem_name[0] == 't':
            return 'magenta'
        elif elem_name[0] == 'f':
            return 'cyan'
