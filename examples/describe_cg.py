#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function, absolute_import, division, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip)


import argparse
from collections import defaultdict, Counter
import logging
import os.path
import itertools as it
import math
import pandas as pd


from logging_exceptions import log_to_exception

import forgi.utilities.commandline_utils as fuc
import forgi.threedee.model.descriptors as ftmd
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.graph_pdb as ftug


log = logging.getLogger(__name__)


def generateParser():
    parser = fuc.get_rna_input_parser("Collect data about a list of rna files and store it as a csv.",
                                      nargs="+", rna_type="any", enable_logging=True)
    parser.add_argument(
        "--csv", type=str, help="Store dataframe under this filename. (Prints to stdout if not given)")
    parser.add_argument("-k", "--keys", type=str, help="Only print the following properties. "
                        "(A comma-seperated list of column headers, e.g. rog_vres")
    parser.add_argument(
        "--angles", type=str, help="Store the angles between the given pairs of elements. Comma-seperated element tuples, seperated by colons. (e.g.: 's0,s1:s1,s2')")
    parser.add_argument("--distances", type=str,
                        help="Store the distances between the given nucleotides. Comma-seperated nucleotide tuples, seperated by colons. (e.g.: '1,20:2,19')")
    parser.add_argument("--per-ml", action="store_true",
                        help="Describe junction segments instead of the whole cg (one entry per segment)")
    parser.add_argument(
        "--mode", type=str, help="For use with --csv. Either 'a' for append or 'o' for overwrite. Default: Abort if file exists.")

    return parser


def describe_rna(cg, file_num, dist_pais, angle_pairs):
    data = {}
    data["nt_length"] = cg.seq_length
    data["num_cg_elems"] = len(cg.defines)
    for letter in "smifth":
        data["num_" + letter] = len([x for x in cg.defines if x[0] == letter])
    multiloops = cg.find_mlonly_multiloops()
    descriptors = []
    junct3 = 0
    junct4 = 0
    reg = 0
    pk = 0
    op = 0
    for ml in multiloops:
        descriptors = cg.describe_multiloop(ml)
        if "regular_multiloop" in descriptors:
            if len(ml) == 3:
                junct3 += 1
            elif len(ml) == 4:
                junct4 += 1
            reg += 1
        if "pseudoknot" in descriptors:
            pk += 1
        if "open" in descriptors:
            op += 1
    data["3-way-junctions"] = junct3
    data["4-way-junctions"] = junct4

    #print (descriptors)
    data["open_mls"] = op
    # print(data["open_mls"][-1])
    data["pseudoknots"] = pk
    data["regular_mls"] = reg
    data["total_mls"] = len(multiloops)
    try:
        data["longest_ml"] = max(len(x) for x in multiloops)
    except ValueError:
        data["longest_ml"] = 0
    try:
        data["rog_fast"] = cg.radius_of_gyration("fast")
    except (ftmc.RnaMissing3dError, AttributeError):
        data["rog_fast"] = float("nan")
        data["rog_vres"] = float("nan")
        data["anisotropy_fast"] = float("nan")
        data["anisotropy_vres"] = float("nan")
        data["asphericity_fast"] = float("nan")
        data["asphericity_vres"] = float("nan")
    else:
        data["rog_vres"] = cg.radius_of_gyration("vres")
        data["anisotropy_fast"] = ftmd.anisotropy(cg.get_ordered_stem_poss())
        data["anisotropy_vres"] = ftmd.anisotropy(
            cg.get_ordered_virtual_residue_poss())
        data["asphericity_fast"] = ftmd.asphericity(cg.get_ordered_stem_poss())
        data["asphericity_vres"] = ftmd.asphericity(
            cg.get_ordered_virtual_residue_poss())
    for from_nt, to_nt in dist_pairs:
        try:
            dist = ftuv.vec_distance(cg.get_virtual_residue(int(from_nt), True),
                                     cg.get_virtual_residue(int(to_nt), True))
        except Exception as e:
            dist = float("nan")
            log.warning("%d%s File %s: Could not calculate distance between "
                        "%d and %d: %s occurred: %s", file_num,
                        {1: "st", 2: "nd", 3: "rd"}.get(
                            file_num % 10 * (file_num % 100 not in [11, 12, 13]), "th"),
                        cg.name, from_nt, to_nt, type(e).__name__, e)
        data["distance_{}_{}".format(from_nt, to_nt)] = dist
    for elem1, elem2 in angle_pairs:
        try:
            angle = ftuv.vec_angle(cg.coords.get_direction(elem1),
                                   cg.coords.get_direction(elem2))
        except Exception as e:
            angle = float("nan")
            log.warning("%d%s File %s: Could not calculate angle between "
                        "%s and %s: %s occurred: %s", file_num,
                        {1: "st", 2: "nd", 3: "rd"}.get(
                            file_num % 10 * (file_num % 100 not in [11, 12, 13]), "th"),
                        cg.name, elem1, elem2, type(e).__name__, e)
        data["angle_{}_{}".format(elem1, elem2)] = angle
    data["missing_residues_5prime"] = (len(cg.seq.with_missing[:1]) - 1)
    data["missing_residues_3prime"] = (
        len(cg.seq.with_missing[cg.seq_length:]) - 1)
    data["missing_residues_middle"] = (
        len(cg.seq.with_missing[1:cg.seq_length]) - len(cg.seq[1:cg.seq_length]))
    data["missing_residues_total"] = (
        len(cg.seq.with_missing[:]) - len(cg.seq[:]))
    fp = len(cg.seq.with_missing[:1]) - 1
    tp = 0
    old_bp = None
    bp = None
    for bp in cg.backbone_breaks_after:
        fp += len(cg.seq.with_missing[bp:bp + 1].split('&')[1]) - 1
        tp += len(cg.seq.with_missing[bp:bp + 1].split('&')[0]) - 1
    tp += len(cg.seq.with_missing[cg.seq_length:]) - 1
    data["missing_residues_5prime_chain"] = (fp)
    data["missing_residues_3prime_chain"] = (tp)
    data["missing_residues_middle_chain"] = (
        data["missing_residues_total"] - fp - tp)
    incomplete_elem_types = Counter(x[0] for x in cg.incomplete_elements)
    data["s_with_missing"] = incomplete_elem_types["s"]
    data["i_with_missing"] = incomplete_elem_types["i"]
    data["m_with_missing"] = incomplete_elem_types["m"]
    data["h_with_missing"] = incomplete_elem_types["h"]
    mp = ""
    if incomplete_elem_types["s"]:
        for elem in cg.incomplete_elements:
            if elem[0] != "s":
                continue
            for i in range(cg.defines[elem][0], cg.defines[elem][1]):
                left_s = cg.seq.with_missing[i:i + 1]
                if len(left_s) > 2:
                    right_s = cg.seq.with_missing[cg.pairing_partner(
                        i + 1):cg.pairing_partner(i)]
                    if len(right_s) > 2:
                        mp += "{}&{};".format(left_s, right_s)
    data["missing_basepairs"] = mp
    return data


def describe_ml_segments(cg):
    data = defaultdict(list)
    loops = cg.find_mlonly_multiloops()
    for loop in it.chain(loops, [[i] for i in cg.iloop_iterator()]):
        print(loop)
        if loop[0][0] == "i":
            description = ["interior_loop"]
        else:
            description = cg.describe_multiloop(loop)
        try:
            j3_roles = cg._assign_loop_roles(loop)
        except ValueError:
            j3_roles = None
        if j3_roles:
            j3_familyFlat = cg._junction_family_westhof1(j3_roles)
            j3_family3D = cg._junction_family_3d(j3_roles)
            j3_familyPerp = cg._junction_family_is_perpenticular(j3_roles)
            j3_Delta = cg.get_length(
                j3_roles["J23"]) - cg.get_length(j3_roles["J31"])
        else:
            j3_family3D = None
            j3_familyFlat = None
            j3_familyPerp = None
            j3_Delta = None

        loop_start = float("inf")
        for segment in loop:
            if cg.define_a(segment)[0] < loop_start:
                loop_start = cg.define_a(segment)[0]
        for segment in loop:
            if segment[0] not in "mi":
                continue
            data["loop_start_after"].append(loop_start)
            data["segment_start_after"].append(cg.define_a(segment)[0])
            data["segment"].append(segment)
            data["junction_length"].append(len(loop))
            data["segment_length"].append(cg.get_length(segment))
            if segment[0] == "i":
                dims = list(sorted(cg.get_bulge_dimensions(segment)))
            else:
                dims = [-1, -1]
            data["iloop_length_1"].append(dims[0])
            data["iloop_length_2"].append(dims[1])
            data["loops_largest_segment_length"].append(
                max(cg.get_length(x) for x in loop))
            data["loops_shortest_segment_length"].append(
                min(cg.get_length(x) for x in loop))
            data["sum_of_loops_segment_lengths"].append(
                sum(cg.get_length(x) for x in loop))
            data["loop_segment_lengths"].append(
                ",".join(map(str, sorted(cg.get_length(x) for x in loop))))

            data["angle_type"].append(
                abs(cg.get_angle_type(segment, allow_broken=True)))
            s1, s2 = cg.connections(segment)

            vec1 = cg.coords.get_direction(s1)
            if cg.get_sides(s1, segment) == (1, 0):
                vec1 = -vec1
            else:
                assert cg.get_sides(s1, segment) == (0, 1)
            vec2 = cg.coords.get_direction(s2)
            if cg.get_sides(s2, segment) == (1, 0):
                vec2 = -vec2
            else:
                assert cg.get_sides(s2, segment) == (0, 1)
            data["angle_between_stems"].append(ftuv.vec_angle(vec1, vec2))
            data["offset1"].append(ftuv.point_line_distance(cg.coords[s1][cg.get_sides(s1, segment)[0]],
                                                            cg.coords[s2][0], cg.coords.get_direction(
                                                                s2)
                                                            ))
            data["offset2"].append(ftuv.point_line_distance(cg.coords[s2][cg.get_sides(s2, segment)[0]],
                                                            cg.coords[s1][0], cg.coords.get_direction(
                                                                s1)
                                                            ))
            closer1, far1 = cg.coords[s1][cg.get_sides(
                s1, segment)[0]], cg.coords[s1][cg.get_sides(s1, segment)[1]]
            closer2, far2 = cg.coords[s2][cg.get_sides(
                s2, segment)[0]], cg.coords[s2][cg.get_sides(s2, segment)[1]]

            data["offset"].append(ftuv.vec_distance(*ftuv.line_segment_distance(closer1, closer1 + (closer1 - far1) * 100000,
                                                                                closer2, closer2 + (closer2 - far2) * 100000)))
            data["junction_va_distance"].append(
                ftug.junction_virtual_atom_distance(cg, segment))
            data["is_external_multiloop"].append("open" in description)
            data["is_pseudoknotted_multiloop"].append(
                "pseudoknot" in description)
            data["is_regular_multiloop"].append(
                "regular_multiloop" in description)
            data["is_interior_loop"].append("interior_loop" in description)
            if j3_roles is not None:
                elem_role, = [x[0]
                              for x in j3_roles.items() if x[1] == segment]
            else:
                elem_role = "?"
            data["j3_role"].append(elem_role)
            data["j3_familyFlat"].append(j3_familyFlat)
            data["j3_family3D"].append(j3_family3D)
            data["j3_familyPerp"].append(j3_familyPerp)
            data["j3_Delta_j23_j31"].append(j3_Delta)
            dssr_stacking = False
            if "dssr_stacks" in cg.infos:
                if segment in cg.infos["dssr_stacks"]:
                    dssr_stacking = True
            data["dssr_stacking"].append(dssr_stacking)

            kh_stem_angle = float("nan")
            if abs(cg.get_angle_type(segment, allow_broken=True)) == 5:
                next_ml = cg.get_next_ml_segment(segment)
                if isinstance(next_ml, str) and next_ml[0] == "m" and abs(cg.get_angle_type(next_ml, allow_broken=True)) == 5:
                    stems1 = cg.edges[segment]
                    stems2 = cg.edges[next_ml]
                    try:
                        s1, s2 = (stems1 | stems2) - (stems1 & stems2)
                    except ValueError:
                        pass
                    else:
                        vec1 = cg.coords.get_direction(s1)
                        vec2 = cg.coords.get_direction(s2)
                        angle = ftuv.vec_angle(vec1, vec2)
                        if angle > math.pi / 2:
                            angle = math.pi - angle
                        kh_stem_angle = angle
            data["kh_stem_angle"].append(kh_stem_angle)
    if data:
        data["pk_number"] = number_by(data, "loop_start_after",
                                      "is_pseudoknotted_multiloop")
        data["loop_number"] = number_by(data, "loop_start_after", None)
        data["reguler_multiloop_number"] = number_by(data, "loop_start_after",
                                                     "is_regular_multiloop")
    return data


def number_by(data, sorting_column="loop_start_after", only_for_col="is_pseudoknotted_multiloop"):
    log.debug(list((key, len(data[key])) for key in data))
    df = pd.DataFrame(data)
    if only_for_col is not None:
        df = df[df[only_for_col] == True]
    sorted_vals = list(sorted(set(df[sorting_column])))
    out_column = []
    for i in range(len(data[sorting_column])):
        if only_for_col is None or data[only_for_col][i]:
            out_column.append(sorted_vals.index(data[sorting_column][i]) + 1)
        else:
            out_column.append(0)
    log.info("number_by column is %s, len(data[%s])=%s)", out_column, sorting_column, len(
        data[sorting_column]))
    return out_column


parser = generateParser()
if __name__ == "__main__":
    args = parser.parse_args()
    cgs, filenames = fuc.cgs_from_args(
        args, "any", enable_logging=True, return_filenames=True)
    data = defaultdict(list)

    if args.distances:
        dist_pairs = str(args.distances).split(str(':'))
        dist_pairs = [x.split(",") for x in dist_pairs]

    else:
        dist_pairs = []
    if args.angles:
        angle_pairs = str(args.angles).split(str(":"))
        angle_pairs = [x.split(",") for x in angle_pairs]
    else:
        angle_pairs = []
    for i, cg in enumerate(cgs):
        file_num = i + 1
        log.info("Describing the %d%s cg %s", file_num, {1: "st", 2: "nd", 3: "rd"}.get(
            file_num % 10 * (file_num % 100 not in [11, 12, 13]), "th"), cg.name)
        try:
            key = {"name": cg.name, "filename": filenames[i]}
            if args.per_ml:
                new_data = describe_ml_segments(cg)
                for i in range(len(new_data["segment"])):
                    for k, v in key.items():
                        data[k].append(v)
                    for k, v in new_data.items():
                        data[k].append(v[i])
            else:
                new_data = describe_rna(cg, file_num, dist_pairs, angle_pairs)
                for k, v in key.items():
                    data[k].append(v)
                for k, v in new_data.items():
                    data[k].append(v)
        except Exception as e:
            with log_to_exception(log, e):
                log.error("Error occurred during describing %d%s cg %s", file_num, {1: "st", 2: "nd", 3: "rd"}.get(
                    file_num % 10 * (file_num % 100 not in [11, 12, 13]), "th"), cg.name)
            raise
    if args.keys:
        allowed_keys = args.keys.split(",") + ["name"]
        for key in list(data.keys()):
            if key not in allowed_keys:
                del data[key]
    df = pd.DataFrame(data)
    df.set_index("name", append=True, inplace=True)
    if args.csv:
        if not args.mode and os.path.isfile(args.csv):
            raise RuntimeError("File {} exists already.".format(args.csv))
        if not args.mode or args.mode == 'o':
            df.to_csv(args.csv)
        elif args.mode == 'a':
            df.to_csv(args.csv, mode='a')
        else:
            raise ValueError("Mode must be one of 'a' and 'o'")
    else:
        pd.set_option('max_columns', 99)
        pd.set_option('max_rows', 200)

        print(df)
