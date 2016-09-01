#!/usr/bin/python

import sys

import itertools as it
import math as m
import numpy as np
import uuid
import collections as col
import warnings

import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.graph_pdb as ftug
#import forgi.threedee.utilities.average_stem_vres_atom_positions as cua
import forgi.utilities.debug as fud
import forgi.threedee.utilities.vector as cuv
import forgi.threedee.utilities.vector as ftuv

import Bio.PDB.Model as bpm
import Bio.PDB.Structure as bps
import Bio.PDB as bp


class PymolPrinter:
    def __init__(self):
        self.display_virtual_residues = False
        self.rainbow = False
        self.basis = None
        self.visualize_three_and_five_prime = True
        self.encompassing_stems = False
        self.state = 2
        self.stem_stem_orientations = None
        self.new_segments = []
        self.segments = []
        self.new_cones = []
        self.cones = []
        self.labels = []
        self.spheres = []
        self.new_spheres = []
        self.boxes = []
        self.virtual_atoms = False
        self.sidechain_atoms = False
        self.override_color = None
        self.element_specific_colors = None
        self.print_text = True
        self.energy_function = None
        self.add_twists = True
        self.add_longrange = False
        self.add_loops = True
        self.chain = None
        self.max_stem_distances = 0
        self.add_letters = False
        self.draw_axes = False
        self.draw_segments = True
        self.pdb_file = None
        self.movie = False
        self.color_modifier=1.0
        self.stem_color = 'green'
        self.multiloop_color = 'red'
        self.prev_obj_name = ''     # The name of the previously created
                                    # object which needs to be hidden
                                    # when creating a movie
        self.only_elements = None

    def get_color_vec(self, color):
        if color == 'green':
            return [0.0, 1.0, 0.0]
        elif color == 'forest':
            return [0.2, 0.6, 0.2]
        elif color == 'blue':
            return [0.0, 0.0, 1.0]
        elif color == 'red':
            return [1.0, 0.0, 0.0]
        elif color == 'orange':
            return [1., 165 / 255., 0.]
        elif color == 'yellow':
            return [1.0, 1.0, 0.0]
        elif color == 'purple':
            return [1.0, 0.0, 1.0]
        elif color == 'white':
            return [1.0, 1.0, 1.0]
        elif color == 'cyan':
            return [0.0, 1.0, 1.0]
        elif color == 'magenta':
            return [249 / 255., 132 / 255., 229 / 255.]
        elif color == 'light gray':
            return [.8, .8, .8]
        elif color == 'dark gray':
            return [.1, .1, .1]
        elif color == 'middle gray':
            return [.6, .6, .6]
        else:
            return [0.0, 0.0, 0.0]

    def add_sphere(self, p, color='green', width=0.2, text="",
                   color_rgb=None):
        if self.override_color is not None:
            color = self.override_color

        if color_rgb is None:
            color_rgb = self.get_color_vec(color)

        self.new_spheres += [(np.array(p), color_rgb, width, text)]

    def transform_spheres(self, translation, rotation):
        for (p, color, width, text) in self.new_spheres:
            p -= translation

            self.spheres += [(p, color, width, text)]

        self.new_spheres = []

    def add_segment(self, p, n, color='green', width=0.2, text="", key=''):

        # exaggerate the length of the stem
        '''
        new_p = p + 3 * cuv.normalize(p - n)
        new_n = n + 3 * cuv.normalize(n - p)

        p = new_p
        n = new_n
        '''
        if self.override_color is not None:
            color = self.override_color

        if not isinstance(color, (list, tuple)):
            color = [str(c * self.color_modifier) for c in self.get_color_vec(color)]        
        else:
            color = [str(c) for c in color[:3]]

        print >>sys.stderr, "color:", color
        #assert(not allclose(p, n))
        self.new_segments += [(np.array(p), np.array(n), color, width, text, key)]

    def add_cone(self, p, n, color='white', width=2.4, text=''):
        if self.override_color is not None:
            color = self.override_color

        cone_extension = 2.
        cyl_vec = cuv.normalize(n-p)
        cyl_len = cuv.magnitude(n-p)

        new_width = width * (cyl_len + cone_extension) / cyl_len

        self.new_cones += [(np.array(p) - cone_extension * cyl_vec, np.array(n), color, width, text)]
        self.new_cones += [(np.array(n) + cone_extension * cyl_vec, np.array(p), color, width, text)]

    def transform_segments(self, translation, rotation):
        for (p, n, color, width, text, key) in self.new_segments:
            p -= translation
            n -= translation
            new_p = np.dot(rotation, p)
            new_n = np.dot(rotation, n)

            self.segments += [(new_p, new_n, color, width, text, key)]

        self.new_segments = []

    def pymol_spheres_string(self):
        self.spheres += self.new_spheres
        s = ''

        for (p, color, width, text) in self.new_spheres:
            color_vec = color
            if np.ndim(p)>0:
                s += "COLOR, %s," % (",  ".join([str(c) for c in color_vec]))
                s += '\n'
                s += "SPHERE, %s, %f," % (", ".join([str(pi) for pi in p]),
                                          width)
                s += '\n'
            else:
                warnings.warn("p is not iterable! It is {} (for '{}'). IGNORING. ".format(p, text))
        return s

    def pymol_axis_string(self):
        w = 0.12  # cylinder width
        l = 10.0  # cylinder length
        h = 3.0  # cone hight
        d = w * 2.618  # cone base diameter
        s = ""

        s += "CYLINDER, 0.0, 0.0, 0.0,   %f, 0.0, 0.0, %f" % (l, w)
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

    def pymol_segments_string(self):
        color = 'green'
        width = 0.2
        s = ''

        self.segments += self.new_segments

        for seg in self.segments:
            (p, n, color, width, text, key) = seg
            if type(color) is not list:
                color_vec = [str(c * self.color_modifier) for c in self.get_color_vec(color)]
            else:
                color_vec = color

            s += " CYLINDER, %f, %f, %f, %f, %f, %f, " % (p[0], p[1], p[2],
                                                          n[0], n[1], n[2])
            s += "%f, %s, %s," % (width, ", ".join(color_vec),
                                  ", ".join(color_vec)) + '\n'
            s += "\n# {}\n".format(key)

        return s

    def pymol_cones_string(self):
        color = 'white'
        width = 0.2
        s = ''
        self.cones += self.new_cones

        for cone in self.cones:
            (p, n, color, width, text) = cone
            color_vec = [str(c) for c in self.get_color_vec(color)]
            s += " CONE, %f, %f, %f, %f, %f, %f, " % (p[0], p[1], p[2],
                                                          n[0], n[1], n[2])
            s += "%f, %s, %s," % (width, ", ".join(color_vec),
                                  ", ".join(color_vec)) + "\n"

        return s

    def pymol_text_string(self):
        counter = 0
        s = ''
        pa_s = 'cmd.set("label_size", 20)\n'
        uids = []

        for (p, n, color, width, text, key) in self.segments:
            if len(text) == 0:
                continue

            # generate a unique identifier for every object so that other
            # scripts can add others that don't clash
            uid = str(uuid.uuid4()).replace('-', 'x')
            uids += [uid]

            s += "cgox_%s = []" % (uid) + '\n'
          
            if np.all(n==p):
                pos = n
                axes = [ [2,0,0], [0,2,0], [0,0,2] ]
            else:
                comp1 = cuv.normalize(n - p)

                ncl = cuv.get_non_colinear_unit_vector(comp1)

                comp2 = cuv.normalize(np.cross(ncl, comp1))
                comp3 = cuv.normalize(np.cross(ncl, comp2))

                pos = (p + n) / 2.0 + 3 * comp2

                axes = [list(comp1 * 2), list(comp2 * 2), list(comp3 * 2)]

            text = "%s: %.1f" % (text, cuv.magnitude(n - p))

            s += "cyl_text(cgox_%s, plain, %s, " % (uid, str(list(pos)))
            s += "\"%s\", 0.20, axes=%s)" % (text, str(axes)) + '\n'
            pa_s += "pa_%s = cmd.pseudoatom(pos=%s," % (uid, str(list(pos)))
            pa_s += "b=1.0, label=\"%s\")\n" % (text)
            counter += 1

        '''
        for (text, pos) in self.labels:
            uid = str(uuid.uuid4()).replace('-', 'x')
            uids += [uid]

            pa_s += "pa_%s = cmd.pseudoatom(pos=%s," % (uid, str(list(pos)))
            pa_s += "b=1.0, label=\"%s\")\n" % (text)
        '''

        s += "cmd.set(\"cgo_line_radius\",0.03)" + '\n'
        for i in range(counter):
            s += "cmd.load_cgo(cgox_%s, " % (uids[i])
            s += "\'cgox%s\')" % (uids[i]) + '\n'
        s += "cmd.zoom(\"all\", 2.0)" + '\n'

        return pa_s

    def pymol_string(self):
        '''
        Output the contents of this structure into a file that can be passed
        in as a pymol script.
        '''

        s = self.pymol_intro_string()

        if self.draw_segments:
            s += self.pymol_segments_string()


        s += self.pymol_spheres_string()

        if self.draw_axes:
            s += self.pymol_axis_string()

        '''
        if self.draw_cones:
            s += self.pymol_cones_string()
        '''

        s += self.pymol_outro_string()

        if self.print_text:
            s += self.pymol_text_string()

        s += self.pymol_box_string()

        return s

    def dump_pdb(self, filename):
        '''
        If the BulgeGraph has a chain created for it, dump that as well.

        @param filename: The filename of the pdb file to which the chain
                         coordinates will be written.
        '''
        if self.chain is None:
            return

        self.chain.child_list.sort()
        mod = bpm.Model(' ')
        s = bps.Structure(' ')

        mod.add(self.chain)
        s.add(mod)

        io = bp.PDBIO()
        io.set_structure(s)
        io.save(filename)

    def dump_pymol_file(self, filename):
        '''
        Output the structure to file.

        @param filename: The location of the output file.
        '''
        # Output the script for showing the coarse-grained elements
        f = open(filename + ".pym", 'w')
        f.write(self.pymol_string())
        f.close()

        # Output the pdb structure
        self.dump_pdb(filename + ".pdb")

        # Output the script file for loading the pdb and coarse grained
        # structure
        f = open(filename + ".pml", 'w')
        f.write("run %s" % (filename + ".pym"))
        f.close()

    def output_pymol_file(self):
        print self.pymol_string()

    def reset(self):
        self.segments = []
        self.new_segments = []
        self.labels = []

    def pymol_intro_string(self):
        self.cgo_uid = str(uuid.uuid4()).replace('-', 'x')
        s = "from pymol.cgo import *" + '\n'
        s += "from pymol import cmd" + '\n'
        s += "from pymol.vfont import plain" + '\n'
        s += "obj%s = [" % (self.cgo_uid) + '\n'
        return s

    def pymol_outro_string(self):
        s = "]" + '\n'

        if self.movie:
            s += "cmd.load_cgo(obj%s, 'ss%s', %d)" % (self.cgo_uid,
                                                      self.cgo_uid,
                                                      self.state) + '\n'
            self.prev_obj_name = self.cgo_uid
            self.state += 1
        else:
            s += "cmd.load_cgo(obj%s, 'ss%s')" % (self.cgo_uid,
                                                  self.cgo_uid) + '\n'

        return s

    def pymol_box_string(self):
        '''
        Pring out the CGO text to describe the boxes.
        '''
        out_str = ''
        for (box, color) in self.boxes:
            uid = str(uuid.uuid4()).replace('-', 'x')
            color_vec = [str(c) for c in self.get_color_vec(color)]
            out_str += 'obj%s = [\n' % (uid)
            out_str += "LINEWIDTH, .8, \n"
            out_str += "BEGIN, LINES, \n"
            out_str += "COLOR, %s," % (",  ".join([str(c) for c in color_vec]))
            out_str += '\n'
            for corner in box:
                out_str += "VERTEX, %f, %f, %f, \n" % (corner[0],
                                                       corner[1],
                                                       corner[2])
            out_str += 'END \n'
            out_str += '] \n'
            out_str += "cmd.load_cgo(obj%s, 'ss%s')\n" % (uid, uid)

        return out_str

    def add_stem_like_core(self, coords, twists, stem_len, key,
                           color='green', width=2.4):
        (p, n) = coords

        self.add_cone(p, n, 'white', width, key)
        self.add_segment(p, n, color, width, key, key=key)
        #self.add_sphere(p, 'light gray', width=2.0 ) 
        #self.add_sphere(n, 'dark gray', width=2.0 ) 

        if self.add_twists:
            mult = 8.
            width = .3
            #twist1o = bg.get_twists(key)[0]
            #twist2o = bg.get_twists(key)[1]
            (twist1o, twist2o) = twists

            self.add_segment(p, p + mult * twist1o, "cyan", width, '', key=key)
            self.add_segment(n, n + mult * twist2o, "magenta", width, '', key=key)
            '''

            twist_rot_mat_l = cuv.rotation_matrix(n - p, -(1.45 / 2.))
            twist_rot_mat_r = cuv.rotation_matrix(n - p, (1.45 / 2.))

            twist1 = np.dot(twist_rot_mat_l, twist1o)
            twist2 = np.dot(twist_rot_mat_l, twist2o)

            twist3 = np.dot(twist_rot_mat_r, twist1o)
            twist4 = np.dot(twist_rot_mat_r, twist2o)



            self.add_segment(p, p + mult * twist1, "white", width, '', key)
            self.add_segment(n, n + mult * twist2, "white", width, '', key)

            self.add_segment(p, p + mult * twist3, "red", width, '', key)
            self.add_segment(n, n + mult * twist4, "red", width, '', key)
            '''

        #stem_len = bg.stem_length(key)

            for i in range(stem_len):
                #(pos, vec) = ftug.virtual_res_3d_pos(bg, key, i)
                res = ftug.virtual_res_3d_pos_core((p, n), twists, i, stem_len)
                (pos, vec_c, vec_l, vec_r) = res
                self.add_segment(pos, pos + mult * vec_c, "orange", width, '', key=key)

                if self.add_letters:
                    self.labels += [('L', list(pos + mult * vec_l))]
                    self.labels += [('R', list(pos + mult * vec_r))]

                '''
                self.add_segment(pos, pos + mult * vec_l, "yellow", width, '', key)
                self.add_segment(pos, pos + mult * vec_r, "purple", width, '', key)
                '''

                if self.display_virtual_residues:
                    self.add_sphere(pos + mult * vec_l, "cyan", 1.)
                    self.add_sphere(pos + mult * vec_r, "magenta", 1.)

        '''
        self.add_sphere(p + mult * twist1, "white", width, key)
        self.add_sphere(n + mult * twist2, "white", width, key)
        '''

    def add_stem_like(self, cg, key, color='green', width=2.4):
        if key in cg.twists:
            return self.add_stem_like_core(cg.coords[key], cg.twists[key],
                                           cg.stem_length(key), key, color, width)
        else:
            return self.add_stem_like_core(cg.coords[key], None,
                                           cg.stem_length(key), key, color, width)

    def draw_bounding_boxes(self, bg, s):
        '''
        Draw bounding boxes for all of the residues encompassed
        by a stem. But only if there is a pdb file handy.

        @param bg: The BulgeGraph
        @param s: The name of the stem
        '''
        if self.pdb_file is None:
            return

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            struct = bp.PDBParser().get_structure('temp', self.pdb_file)
        chain = list(struct.get_chains())[0]

        for i in range(bg.stem_length(s)):
            (origin, bases, bb) = ftug.bounding_boxes(bg, chain, s, i)
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
                    new_corners += [origin + cuv.change_basis(np.array(corner),
                                    cuv.standard_basis, bases[k])]
                corners = np.array(new_corners)

                if k == 0:
                    self.boxes += [(corners, 'yellow')]
                    self.add_sphere(corners[0], 'yellow', 0.4, '',
                                    [238 / 255., 221 / 255., 130 / 255.])
                    self.add_sphere(corners[7], 'yellow', 0.4, '',
                                    [184 / 255., 134 / 255., 11 / 255.])
                else:
                    self.add_sphere(corners[0], 'purple', 0.4, '',
                                    [238 / 255., 130 / 255., 238 / 255.])
                    self.add_sphere(corners[7], 'purple', 0.4, '',
                                    [208 / 255., 32 / 255., 144 / 255.])
                    self.boxes += [(corners, 'purple')]

    def add_encompassing_cylinders(self, cg, radius=7.):
        cylinders_to_stems = ftug.get_encompassing_cylinders(cg, radius)

        for stems in cylinders_to_stems.values():
            print "stems:", stems

            points = []
            for s in stems:
                points += [cg.coords[s][0], cg.coords[s][1]]
            
            # create the linear regression
            data = np.array(points)
            datamean = data.mean(axis=0)

    
            uu, dd, vv = np.linalg.svd(data - datamean)

            furthest = max([ftuv.magnitude(d) for d in (data - datamean) ])

            start_point = -furthest * vv[0] + datamean
            end_point = furthest * vv[0] + datamean

            self.add_segment(start_point, end_point, 'white', width=4, text='', key=key)


        print >>sys.stderr, "YOOOOOOOOOOOOOOOOOOOOOOO"

    def get_element_color(self, elem_name):
        '''
        Get the color for this element. The color is determined by the name
        of the element.

        @param elem_name: The name of the element.
        @return: A string with a color name
        '''
        if self.element_specific_colors is not None:
            if elem_name in self.element_specific_colors:
                return self.element_specific_colors[elem_name]

        print >>sys.stderr, "self.stem_color:", self.stem_color
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

    def add_dashed(self, point1, point2, width=0.3):
        '''
        Add a dashed line from point1 to point2.
        '''
        dash_length = 0.6
        gap_length = dash_length * 2
        direction = ftuv.normalize(point2 - point1)

        num_dashes = ftuv.magnitude(point2 - point1) / (dash_length + gap_length)
        key=None

        for i in range(int(num_dashes)):
            self.add_segment(point1 + i * (dash_length + gap_length) * direction, 
                             point1 + (i * (dash_length + gap_length) + dash_length) * direction, "purple",
                             width, "", key=key)


    def coordinates_to_pymol(self, cg):
        loops = list(cg.hloop_iterator())

        for key in cg.coords.keys():
            if self.only_elements is not None:
                if key not in self.only_elements:
                    continue

            (p, n) = cg.coords[key]
            color = self.get_element_color(key)

            if key[0] == 's':
                self.add_stem_like(cg, key, color=color)
                self.draw_bounding_boxes(cg, key)
            else:
                #self.add_sphere(p+(n-p)*0.2, 'light gray', width=1.5 ) 
                #self.add_sphere(n+(p-n)*0.2, 'dark gray', width=1.5 )
                if key[0] == 'h':
                    if self.add_loops:
                        if key in loops:
                            self.add_segment(p, n, color, 1.0,
                                             key + " " + str(cg.get_length(key)),
                                            key=key)
                elif key[0] == 'm':
                    twists = cg.get_twists(key)

                    # check if the multiloop is longer than one. If it's not, then
                    # it has an empty define and we its length will be 1
                    if len(cg.defines[key]) == 0:
                        self.add_segment(p, n, color, 1.0,
                                         key + " 1", key=key)
                    else:
                        self.add_segment(p, n, color, 1.0,
                                         key + " " +
                                         str(cg.defines[key][1] -
                                         cg.defines[key][0] + 1), key=key)

                    self.add_segment(p, p+ 7 * twists[0], 'light gray', 0.3, '', key=key)
                    self.add_segment(n, n+ 7 * twists[1], 'light gray', 0.3, '', key=key)

                    x = (p + n) / 2
                    t = ftuv.normalize((twists[0] + twists[1]) / 2.)
                    self.add_segment(x, x + 7 * t, 'middle gray', 0.3, key=key)
                elif key[0] == 'f':
                    if self.visualize_three_and_five_prime:
                        self.add_segment(p, n, color, 1.0,
                                         key + " " +
                                         str(cg.defines[key][1] -
                                         cg.defines[key][0] + 1) + "", key=key)

                elif key[0] == 't':
                    if self.visualize_three_and_five_prime:
                        self.add_segment(p, n, color, 1.0,
                                         key + " " +
                                         str(cg.defines[key][1] -
                                         cg.defines[key][0]) + "", key=key)
                else:
                    #self.add_stem_like(cg, key, "yellow", 1.0)
                    self.add_segment(p, n, color, 1.0, key, key=key)

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

                        self.add_dashed(p[0], p[1])

                        '''
                        self.add_segment(point1, point2, "purple",
                                         0.3, key1 + " " + key2, key=key)
                        
                        '''
                    except:
                        continue

        if self.encompassing_stems:
            self.add_encompassing_cylinders(cg, 7.)

        if self.max_stem_distances > 0:
            for (s1, s2) in it.permutations(cg.stem_iterator(), r=2):
                (i1, i2) = cuv.line_segment_distance(cg.coords[s1][0],
                                                     cg.coords[s1][1],
                                                     cg.coords[s2][0],
                                                     cg.coords[s2][1])
                if cuv.magnitude(i2 - i1) < self.max_stem_distances:
                    #self.add_segment(i1, i2, 'cyan', 0.3, s1 + " " + s2, key=key)
                    self.add_segment(i1, i2, 'cyan', 0.3, key=key)

        if self.virtual_atoms or self.sidechain_atoms:
            cg.add_all_virtual_residues()
            va = ftug.virtual_atoms(cg, sidechain=self.sidechain_atoms)

            atom_width = 0.5
            for i,r in enumerate(sorted(va.keys())):
                for a in va[r].keys():
                    if self.rainbow:
                        import matplotlib
                        matplotlib.use('Agg')
                        import matplotlib.pyplot as plt
                        cmap = plt.get_cmap('gist_rainbow')
                        self.add_sphere(va[r][a], 
                                        color_rgb = cmap(i / float(len(va.keys()))), 
                                        width=atom_width)
                    else:
                        d = cg.get_node_from_residue_num(r)
                        if d[0] == 's':
                            if a in ftup.nonsidechain_atoms:
                                self.add_sphere(va[r][a], self.stem_color, width=atom_width)
                            else:
                                self.add_sphere(va[r][a], 'forest', width=atom_width)
                        elif d[0] == 'i':
                            self.add_sphere(va[r][a], 'yellow', width=atom_width)
                        elif d[0] == 'm':
                            self.add_sphere(va[r][a], self.multiloop_color, width=atom_width)
                        elif d[0] == 'h':
                            self.add_sphere(va[r][a], 'blue', width=atom_width)

        if self.basis:
            for d in cg.defines.keys():
                origin, basis = ftug.element_coord_system(cg, d)

                self.add_segment(origin, origin + 7. * basis[1], 'purple', 2., key=key)

        print >>sys.stderr, "energy_function:", self.energy_function
        # print the contributions of the energy function, if one is specified
        if self.energy_function is not None:
            print >>sys.stderr, "key"
            sum_energy = 0.

            e_func = self.energy_function
            e_func_iter = e_func.interaction_energy_iter(cg, background=False)
            int_energies = list(e_func_iter)
            max_energy = max(int_energies, key=lambda x: x[1])
            print >>sys.stderr, "max_energy:", max_energy

            for (interaction, energy) in int_energies:
                (p, n) = (cg.get_point(interaction[0]),
                          cg.get_point(interaction[1]))
                scaled_energy = - max_energy[1] + energy

                self.add_segment(p, n, 'purple', 3 * np.exp(scaled_energy), key=key)

                sum_energy += energy

        if self.stem_stem_orientations is not None:
            for (s1, s2) in it.permutations(cg.stem_iterator(), 2):
                '''
                if cg.are_adjacent_stems(s1, s2):
                    continue
                '''

                if s1 != 's65':
                    if s2 != 's65':
                        continue

                s1_vec = cg.coords[s1][1] - cg.coords[s1][0]
                s2_vec = cg.coords[s2][1] - cg.coords[s2][0]
                (i1, i2) = cuv.line_segment_distance(cg.coords[s1][0],
                                                     cg.coords[s1][1],
                                                     cg.coords[s2][0],
                                                     cg.coords[s2][1])
                i_vec = i2 - i1

                #i_rej will be orthogonal to s1_vec in the direction
                #of s2
                i_rej = cuv.vector_rejection(i_vec, s1_vec)

                #plane_vec will be orthogonal to s1_vec and to the direction
                # of s2
                plane_vec = np.cross(i_rej, s1_vec)

                # s2_proj is in the intersection plane
                s2_proj_in = cuv.vector_rejection(s2_vec, plane_vec)
                # s2 proj_out is out of the intersection plane
                #s2_proj_out = cuv.vector_rejection(s2_vec, i_rej)

                start_point = cg.coords[s1][0] + 5 * cg.twists[s1][0]
                ortho_offset = cuv.magnitude(i_rej)
                dist = cuv.magnitude(i_vec) + 0.0001

                lateral_offset = m.sqrt(dist ** 2 - ortho_offset ** 2)

                if lateral_offset > 10:
                    continue

                '''
                #self.add_segment(start_point,
                                  start_point + 10 * cuv.normalize(s2_vec),
                                  'white', 0.5)
                #self.add_segment(start_point,
                                  start_point + 5 * cuv.normalize(plane_vec),
                                  'magenta', 0.5)
                #self.add_segment(start_point,
                                  start_point + 5 * cuv.normalize(i_vec),
                                  'cyan', 0.5)
                #self.add_segment(i1, i1 + i_rej,  'cyan', 0.5)
                '''
                self.add_segment(start_point,
                                 start_point + 7 * cuv.normalize(s2_proj_in),
                                 'white', 1.5, key=key)
                '''
                #self.add_segment(start_point,
                                  start_point + 7 * cuv.normalize(s2_proj_out),
                                  'blue', 0.5)
                '''

    def chain_to_pymol(self, chain):
        '''
        Add a Bio.PDB.Chain to the structure, so that it can later be printed.
        '''
        self.chain = chain

    def load_flex_stats(self, flex_file):
        f = open(flex_file, 'r')

        d = col.defaultdict(lambda: col.defaultdict(float))

        for line in f.readlines():
            parts = line.strip().split(' ')

            d[int(parts[0])][int(parts[1])] = float(parts[2])

        return d

    def flex_to_pymol(self, cg, flex_file):
        flex_stats = self.load_flex_stats(flex_file)

        for key in cg.defines.keys():
            if key[0] != 's':
                if key in cg.coords.keys():
                    coords = cg.coords[key]
                    p = (coords[1] + coords[0]) / 2.

                    bd = cg.defines[key]

                    if len(bd) == 2:
                        #out_str += "0 %d" % (abs(bd[1] - bd[0]) + 1)
                        dims = (0, abs(bd[1] - bd[0]) + 1)
                    else:
                        dims = (abs(bd[1] - bd[0]) + 1, abs(bd[2] - bd[3]) + 1)
                        #out_str += "%d %d" % ( min(dims), max(dims))

                    flex = flex_stats[min(dims)][max(dims)] * 10.

                    if key[0] == 'm':
                        self.add_sphere(p, "red", flex, key)
                    elif key[0] == 'i':
                        self.add_sphere(p, "yellow", flex, key)

    def centers_to_pymol(self, cg):

        for key in cg.defines.keys():
            if key in cg.coords.keys():
                coords = cg.coords[key]
                p = (coords[1] + coords[0]) / 2.

                if key[0] == 's':
                    self.add_sphere(p, 'green', 3, key)
                else:
                    if len(cg.edges[key]) == 1:
                        self.add_sphere(p, 'blue', 1.5, key)
                if len(cg.edges[key]) == 2:
                    if key[0] == 'm':
                        self.add_sphere(p, "red", 1.5, key)
                    else:
                        self.add_sphere(p, "yellow", 1.5, key)

    """def stem_atoms(self, coords, twists, stem_len, side=0):
        '''
        Add the locations of the virtual atoms as spheres.

        @param coords: The start and end coordinates of the stem.
        @param twists: The two twists of the stem.
        @param stem_len: The length of the stem.
        '''
        prev_p = [None, None]
        first_p = [None, None]
        last_o3 = [None, None]
        first_o3 = [None, None]

        colors = ['yellow', 'purple']

        for i in range(stem_len):
            vbasis = ftug.virtual_res_basis_core(coords, twists, i, stem_len)
            vpos = ftug.virtual_res_3d_pos_core(coords, twists, i, stem_len)

            # iterate once for each strand
            j = side
            # just use A for now
            for a in cua.avg_stem_vres_atom_coords[j]['A'].items():
                c = a[1]
                new_coords = np.dot(vbasis.transpose(), c) + vpos[0]
                #self.add_sphere(new_coords, colors[j], 0.3)

                if a[0] == 'P' and i == 0:
                    first_p[j] = new_coords
                if a[0] == 'P':
                    if prev_p[j] is not None:
                        self.add_segment(prev_p[j], new_coords,
                                         colors[j], 0.7)
                    prev_p[j] = new_coords
                if a[0] == 'O3*' and i == 0:
                    first_o3[j] = new_coords
                if a[0] == 'O3*':
                    last_o3[j] = new_coords

        #self.add_segment(prev_p[0], last_o3[0], colors[0], 0.7)
    """
