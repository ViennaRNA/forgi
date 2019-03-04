#!usr/bin/env python
from __future__ import division

import logging
import argparse
import copy
import csv
from collections import defaultdict
import pandas

import forgi.graph.bulge_graph as fgb
import forgi.utilities.commandline_utils as fuc
import forgi.threedee.utilities.vector as ftuv
from forgi.utilities.exceptions import GraphIntegrityError
from forgi.utilities.exceptions import GraphConstructionError
from forgi.utilities.numbered_dotbracket import NumberedDotbracket

log = logging.getLogger(__name__)

"""
Classify all pseudoknots found in the input 3D structure based on their shapes
according to Reidys et al. concept into H-type, kissing hairpin,...

Additionally, calculate the angles between stems in the poseudoknots.
Annotation Tools: MC-mcannotate or DSSR (enable stacking analysis)

See:
Reidys, C.M. et al., 2011. Topology and prediction of RNA pseudoknots.
Bioinformatics, 27(8), pp.1076-1085.

Input:
python pseudoknot_analyzer.py filename.cg --minlength #

:param minlength: Minimum length of each stem, default=2
"""

PK_CHECK = ["[","]","{","}","<",">"]
CLOSE_TO_OPEN = {
        ")":"(",
        "]":"[",
        "}":"{",
        ">":"<"
        }

def identification_pseudoknot(rna, without_pk, with_pk, unfold, minlength):
    """
    Return isolated pseudoknots from the dotbracket string

    :param rna: A forgi CoarseGrainRNA object
    :param without_pk: Counter for structures without a pseudoknots
    :param with_pk: Counter for structures with a pseudoknots
    :param unfold: Counter for unfolded Structures
    """
    if isinstance(rna, fgb.BulgeGraph):
        residue_numbers = copy.copy(rna.seq._seqids)
        dotbracket = NumberedDotbracket(rna.to_dotbracket_string(), residue_numbers)
    else:
        dotbracket = rna

    #variables
    modified_dotbracket = NumberedDotbracket()
    pseudoknots_dictionary = {}
    counter_pseudoknots = 0
    pseudoknots = []
    pseudoknot = NumberedDotbracket()
    #general, round, square, curve, arrow
    pseudoknot_open = [0,0,0,0,0]
    pseudoknot_check = False

    #remove all dots:
    dotbracket = dotbracket.without_unpaired()
    log.info("Bracket-structure: {}".format(dotbracket))

    if not dotbracket:
        unfold += 1
        return  without_pk, with_pk, unfold, pseudoknots
    else:
        pass

    dotbracket = dotbracket.without_short_helices(minlength)
    dotbracket = dotbracket.without_unpaired()

    if not dotbracket:
        unfold += 1
        return  without_pk, with_pk, unfold, pseudoknots
    else:
        pass

    #remove consecutive basepairs = "condense"
    dotbracket  = dotbracket.condensed()
    log.info("Condensed-bracket-structure: {}".format(dotbracket))

    if set(dotbracket)-set("()[]{}<>"):
        raise ValueError("Pseudoknots with more than 4 types of "
                         "brackets are not supported.")

    #discard every structure without pseudoknot
    if not any (pk in dotbracket for pk in PK_CHECK):
        without_pk += 1
        print("Without Psedudoknot: {}".format(without_pk))
    #for structures with pseudoknots
    else:
        #for number, residue in residue_numbers.items():
        with_pk += 1
        for bracket in dotbracket:
            log.info("actual bracket {}".format(bracket))
            #CLOSE BRACKETS
            if bracket in CLOSE_TO_OPEN:
                log.info("CLOSE pk-check {}".format(pseudoknot_check))
                #close pseudoknot
                if pseudoknot_check == True:
                    openbracket = CLOSE_TO_OPEN[bracket]
                    log.debug("bracket:{}; pseudoknot[-1]:{}; openbracket: {} "
                    .format(bracket, pseudoknot[-1],openbracket))

                    #update open brackets
                    pseudoknot_open[0] -= 1
                    if bracket == PK_CHECK[1]:
                        pseudoknot_open[2] -= 1
                    elif bracket == PK_CHECK[3]:
                        pseudoknot_open[3] -= 1
                    elif bracket == PK_CHECK[5]:
                        pseudoknot_open[4]-= 1
                    else:
                        pseudoknot_open[1]-= 1

                    #pseudoknot uptdate
                    if openbracket == pseudoknot[-1]:
                        pseudoknot = pseudoknot[:-1]
                        log.info("check: {}; pseudoknot yet: {}; open: {}"
                        .format(pseudoknot_check, pseudoknot, pseudoknot_open))
                    else:
                        pseudoknot = pseudoknot + bracket
                        #check if the pseudoknot is finished
                        if not any(pseudoknot_open):
                            pseudoknots.append(pseudoknot)
                            counter_pseudoknots += 1
                            pseudoknots_dictionary.update({counter_pseudoknots:pseudoknot})
                            pseudoknot_check = False
                            pseudoknot = NumberedDotbracket()
                            pseudoknot_open = [0,0,0,0,0]
                            log.info("FINISHED PSEUDOKNOT: check: {}; pseudoknot yet: {}; pseudoknots: {}"
                            .format(pseudoknot_check, pseudoknot, pseudoknots))
                        #check if bracket from modified_dotbracket is needed
                        elif pseudoknot_open[0]==0 or any(v < 0 for v in pseudoknot_open):
                            pseudoknot = modified_dotbracket[-1] + pseudoknot
                            modified_dotbracket = modified_dotbracket[:-1]
                            pseudoknot_open[0] += 1
                            pseudoknot_open[1] += 1
                            log.info("modified_dotbracket:{}; pseudoknot yet: {}; open: {} "
                            .format(modified_dotbracket, pseudoknot, pseudoknot_open))
                        #continue
                        else:
                            log.info("check: {}; pseudoknot yet: {}; open: {}"
                            .format(pseudoknot_check, pseudoknot, pseudoknot_open))

                #normal close
                elif modified_dotbracket[-1] == CLOSE_TO_OPEN[bracket]:
                    modified_dotbracket = modified_dotbracket[:-1]
                    log.info("skip last modified_dotbracket {}".format(modified_dotbracket))
                else:
                    log.info("Something strange happen - close function")
                    assert False

            #OPEN BRACKETS
            #if a pesoudoknot starts [ (can only start with "[")
            elif (bracket == PK_CHECK[0]) and pseudoknot_check == False:
                #position before can only be a (, all ) should be already skipped
                pseudoknot = modified_dotbracket[-1] + bracket
                modified_dotbracket = modified_dotbracket[:-1]
                pseudoknot_check = True
                pseudoknot_open[0] += 2 #([
                pseudoknot_open[1] += 1 #(
                pseudoknot_open[2] += 1 #[
                log.info("PK-START: modified_dotbracket: {}; check: {}; pseudoknot jet: {}, open {}"
                            .format(modified_dotbracket, pseudoknot_check, pseudoknot, pseudoknot_open))
            else:
                #if PK is True
                if pseudoknot_check == True:
                    pseudoknot += bracket
                    pseudoknot_open[0] += 1
                    if bracket == PK_CHECK[0]:
                        pseudoknot_open[2] += 1
                    elif bracket == PK_CHECK[2]:
                        pseudoknot_open[3] += 1
                    elif bracket == PK_CHECK[4]:
                        pseudoknot_open[4] += 1
                    else:
                        pseudoknot_open[1] += 1
                    log.info("check: {}; pseudoknot yet: {}, open: {}"
                    .format(pseudoknot_check, pseudoknot, pseudoknot_open))
                #if PK is False
                else:
                    modified_dotbracket +=bracket
                    log.info("open: {}; append modified_dotbracke: {}"
                    .format(bracket,modified_dotbracket))

    if not modified_dotbracket:
        log.info("modified_dotbracket is empty - WORKS CORRECT HERE")
    else:
        log.info("Something went wrong due skipping () pairs, rest: {}"
        .format(modified_dotbracket))
        assert False

    return without_pk, with_pk, unfold, pseudoknots


def classify_pseudoknot_genus1(pseudoknot):
    """
    Return the classes of a pseudoknot with genus1

    :param pseudoknot: Current genus1-"other" pseudoknot
    """
    if pseudoknot == "([)]":
        return "g1_H_type"
    elif pseudoknot == "([)(])":
        return "g1_Kissinghairpin"
    elif pseudoknot == "([{)]}":
        return "g1_L"
    elif pseudoknot == "([{)(]})":
        return "g1_M"
    else:
        log.debug("Other pseudoknot: %s", str(pseudoknot))
        return "other"

def update_pk_classes(classes1, classes2):
    """
    Update the count of the current pseudoknot class

    :param classes1: Dictionary with pseudoknot classes genus1
    :param classes2: Dictionary with pseudoknot classes genus2
    """
    for k,v in classes2.items():
        classes1[k]+=v

def classify_pseudoknots(pseudoknots):
    """
    Return the class of a pseudoknot with and counts them

    :param pseudoknots: A list of pseudoknots
    """
    pk_classes = {
        "g1_H_type":[],
        "g1_Kissinghairpin":[],
        "g1_L":[],
        "g1_M":[],
        "g2_H_type":[],
        "g2_Kissinghairpin":[],
        "g2_L":[],
        "g2_M":[],
    }
    other = []
    for pseudoknot in pseudoknots:
        log.debug("Classifying %s for genus 1", pseudoknot)
        pk_type = classify_pseudoknot_genus1(pseudoknot)
        #if the pseudoknot is no other - append it
        if pk_type != "other":
            pk_classes[pk_type].append(pseudoknot)
        #if the pseudoknot is "other" clean it again
        else:
            log.debug("Other %s", pseudoknot)
            _,_,_,cleaned_pks = identification_pseudoknot(pseudoknot,0,0,0,0)
            #if cleand version is different classify again
            if cleaned_pks != [pseudoknot]:
                classes_cleaned, cleaned_other = classify_pseudoknots(cleaned_pks)
                other+=cleaned_other
                update_pk_classes(pk_classes, classes_cleaned)
            #otherwise same with genus 2
            else:
                cleaned_pk, = cleaned_pks #Almost same as cleaned_pk = cleaned_pks
                classes_g2, remaining_pks = classify_pseudoknot_genus2(cleaned_pk)
                other+=remaining_pks
                update_pk_classes(pk_classes, classes_g2)
    return pk_classes, other

def classify_pseudoknot_genus2(pseudoknot):
    """
    Return the classes of a pseudoknot with genus2 and count the occurences

    :param pseudoknot: Current genus1-"other" pseudoknot
    """
    log.debug("Classifying %s for genus 2", pseudoknot)
    pk_classes = {
        "g2_H_type":[],
        "g2_Kissinghairpin":[],
        "g2_L":[],
        "g2_M":[],
    }
    while True:
        old_pseudoknot = pseudoknot
        if "([)]" in pseudoknot:
            occurrences, pseudoknot = pseudoknot.without_substr("([)]")
            pk_classes["g2_H_type"]+=occurrences
        if "([)(])" in pseudoknot:
            occurrences, pseudoknot = pseudoknot.without_substr("([)(])")
            pk_classes["g2_Kissinghairpin"]+=occurrences
        if "([{)]}" in pseudoknot:
            occurrences, pseudoknot = pseudoknot.without_substr("([{)]}")
            pk_classes["g2_L"]+=occurrences
        if pseudoknot == "([{)(]})":
            occurrences, pseudoknot = pseudoknot.without_substr("([{)(]})")
            pk_classes["g2_M"]+=occurrences
        if pseudoknot == old_pseudoknot:
            break
    if pseudoknot:
        _,_,_,cleaned_pks = identification_pseudoknot(pseudoknot,0,0,0,0)
        remaining_pks = []
        if cleaned_pks != [pseudoknot]:
            for pk in cleaned_pks:
                pk_classes_2, remaining = classify_pseudoknot_genus2(pk)
                update_pk_classes(pk_classes, pk_classes_2)
                remaining_pks+=remaining
        else:
            remaining_pks = [pseudoknot]
    else:
        remaining_pks = []
    return pk_classes, remaining_pks


def stem_parameters(stem, rna, side):
    """
    :param stem: E.g. 's1'
    :param rna: A forgi CoarseGrainRNA object
    :param side: 0 or 1. See explaination above.
    """
    position = rna.coords[stem][side]
    direction = rna.coords[stem][not side]-position
    return position, direction

def stem_after_next_ml(rna, pos, before):
    ml_found = False
    for elem in rna.iter_elements_along_backbone(pos):
        if elem[0]=="m":
            ml_found = True
        elif ml_found:
            assert elem[0]=="s"
            return elem
        elif elem == before:
            return None
    return None

def extend_pk_description(dataset, filename, pk_type, rna, pk, pk_number):
    """
    Return a extended descripiton of current pseudoknot in the current files
    e.g. angles between stems

    :param dataset:  Current dataset that will be updated
    :param filename: Filename of the current structure
    :parma pk_type: Class of the pseudoknot
    :param rna: A forgi CoarseGrainRNA object
    :param pk: Structure of the pseudoknot, a NumberedDotbracket object,
               in a condensed (shadow-like) representation.
               This representation always contains the most 5' basepair.
    :param pk_number: consecutive number of the pseudoknot
    """
    domains = rna.get_domains()
    helices = domains["rods"] # A list of elements, e.g. ["s0", "i0", "s1"]
    log.debug("Helices: %s", helices)
    #rna.log(logging.WARNING)
    stems_5p = []
    stems_3p = []

    nums = []
    log.debug("pk Residue numbers %s", pk.residue_numbers)
    log.debug("pk helix ends %s", pk.helix_ends)

    for i, resnum in enumerate(pk.residue_numbers):
        num = rna.seq.to_integer(resnum)
        nums.append(num)
        element_5p = rna.get_node_from_residue_num(num)
        stems_5p.append(element_5p)

        num2 = rna.seq.to_integer(pk.helix_ends[i])
        log.debug("num %s nums2 %s", num, num2)
        element_3p =rna.get_node_from_residue_num(num2)
        stems_3p.append(element_3p)
    log.debug("nums %s", nums)
    for i, stem1_5p in enumerate(stems_5p):
        dataset["Filename"].append(filename)
        dataset["rnaname"] = rna.name
        dataset["pk_type"].append(pk_type)
        dataset["pk_id"].append(pk_number)
        dataset["angle_nr"].append(i)
        if pk_type == "other":
            dataset["pk_structure"].append(str(pk))
        else:
            dataset["pk_structure"].append("")
        #is this the first occurrence of stem in stems?
        if stems_5p.index(stem1_5p)==i:
            #first occurrence. Strand 0, look at 3' end of helix
            stem1 = stems_3p[i]
            strand = 0
        else:
            assert i>stems_5p.index(stem1_5p)
            stem1 = stem1_5p
            strand = 1
        try:
            stem2_5p = stems_5p[i+1]
        except IndexError:
            stem2_5p = stems_5p[0]
            outside_pk = True
        else:
            outside_pk = False
        if outside_pk or stems_5p.index(stem2_5p)==i+1:
            #first occurrence
            stem2 = stem2_5p
            strand2 = 0
        else:
            strand2 = 1
            if outside_pk:
                stem2 = stems_3p[0]
            else:
                stem2 = stems_3p[i+1]
        log.debug("Stem 5' %s, 3' %s, stem1 %s stem2 %s", stems_5p, stems_3p, stem1, stem2)
        # enable stacking analysis via DSSR
        # differentiate between stacking (True), no stacking (False) and brakes
        # within/aorund the pseudoknot (-1) incl. 'virtual' angles e.g. H-Type angle_type3
        ml_stack=[]
        if rna.dssr:
            nc_bps = list(rna.dssr.noncanonical_pairs())
            nc_dict = defaultdict(list)
            for nt1, nt2, typ in nc_bps:
                nc_dict[nt1].append((nt2, typ))
                nc_dict[nt2].append((nt1, typ))
            stacking_loops = rna.dssr.stacking_loops()
            start_found = 0
            connection = []
            stacking = None
            branch = None
            log.debug("Checking %s and %s for stacking, strand %s", stem1, stem2, strand)
            for elem in rna.iter_elements_along_backbone(): #walk along the backbone
                if start_found == strand+1:
                    if branch:
                        log.debug("in branch: elem %s, branch %s, stacking %s", elem, branch, stacking)
                        if elem == branch:
                            log.debug("End branch at %s", elem)
                            branch = None
                            log.debug("Branch end")
                        continue
                    if elem[0] != "s":
                        connection.append(elem)
                        if rna.defines[elem] and rna.defines[elem][-1] in rna.backbone_breaks_after:
                            stacking = -1
                        if elem not in stacking_loops and stacking != -1:
                            stacking = False
                    elif elem == stem2:
                        if stacking is None:
                            stacking = True
                        log.debug("Found second stem, elem %s, stacking %s", elem, stacking)
                        break
                    elif elem[0] == "s" and connection:
                        branch = elem
                        if rna.defines[elem][-1] in rna.backbone_breaks_after:
                            stacking = -1
                    log.debug("elem %s, stacking %s, branch %s", elem, stacking, branch)
                elif elem == stem1:
                    start_found += 1
                    if rna.defines[elem][strand*2+1] in rna.backbone_breaks_after:
                        stacking = -1
                    log.debug("First stem, elem %s, stacking %s", elem, stacking)
            else:
                log.debug("End iteration, stacking->-1")
                stacking = -1
            log.debug("Finally, stacking = %s", stacking)
            # more detailed stacking (including backbone brackes within and around the pseudoknot)
            dataset["this_loop_stacking_dssr"].append(stacking)
            dataset["connecting_loops"].append(",".join(connection))

            # more genereal stacking information
            connecting_loops = rna.edges[stem1]&rna.edges[stem2]
            for loop in connecting_loops:
                if loop in stacking_loops:
                    ml_stack.append(loop)
            stacks = rna.dssr.coaxial_stacks()
            log.info("Stacks: %s", stacks)
            for stack in stacks:
                if stem1 in stack and stem2 in stack:
                    # the two stems stack, but we do not specify along which
                    # multiloop segment they stack.
                    dataset["is_stacking_dssr"].append(True)
                    break
            else:
                dataset["is_stacking_dssr"].append(False)

            # Does the connection form base-triples with the stem?
            stem1_triples=0
            stem2_triples=0
            aminors1 = 0
            aminors2 = 0
            aminors = list(rna.dssr.aminor_interactions())
            for elem in connection:
                for nt in rna.define_residue_num_iterator(elem,seq_ids=True):
                    if (nt, stem1) in aminors:
                        aminors1+=1
                        log.debug("AMinor %s (%s), %s", nt, elem, stem1)
                    elif (nt, stem2) in aminors:
                        aminors2+=1
                        log.debug("AMinor %s (%s), %s", nt, elem, stem2)
                    else:
                        for partner, typ in nc_dict[nt]:
                            if rna.get_elem(partner)==stem1:
                                log.debug("base_triple %s, %s: %s-%s (%s)", elem, stem1, nt,partner,typ)
                                stem1_triples+=1
                            elif rna.get_elem(partner)==stem2:
                                log.debug("base_triple %s, %s: %s-%s (%s)", elem, stem2, nt,partner,typ)
                                stem2_triples+=1
            log.debug("%s has a length of %s and %s triples", stem1, rna.stem_length(stem1),stem1_triples)
            log.debug("%s has a length of %s and %s triples", stem2, rna.stem_length(stem2),stem2_triples)
            dataset["stem1_basetripleperc_dssr"].append(stem1_triples/rna.stem_length(stem1))
            dataset["stem2_basetripleperc_dssr"].append(stem2_triples/rna.stem_length(stem2))
            dataset["stem1_aminorperc_dssr"].append(aminors1/rna.stem_length(stem1))
            dataset["stem2_aminorperc_dssr"].append(aminors2/rna.stem_length(stem2))

        else:
            dataset["is_stacking_dssr"].append(float("nan"))
            dataset["this_loop_stacking_dssr"].append(float("nan"))
            dataset["connecting_loops"].append("")
            dataset["stem1_basetripleperc_dssr"].append(float("nan"))
            dataset["stem2_basetripleperc_dssr"].append(float("nan"))
            dataset["stem1_aminorperc_dssr"].append(float("nan"))
            dataset["stem2_aminorperc_dssr"].append(float("nan"))

        dataset["stacking_loops"].append(",".join(ml_stack))

        pos1, dir1 = stem_parameters(stem1, rna, not strand)
        pos2, dir2 = stem_parameters(stem2, rna, strand2)
        dataset["stem1"].append(stem1)
        dataset["stem2"].append(stem2)

        dataset["angle_between_stems"].append(ftuv.vec_angle(dir1, dir2))
        dataset["distance_between"].append(ftuv.vec_distance(pos1, pos2))

        next_stem = None
        if not outside_pk:
            next_stem = stem_after_next_ml(rna, nums[i], before=stem2)
            if next_stem==stem2:
                next_stem  = None
        if next_stem:
            posN, dirN = stem_parameters(next_stem, rna, 0)
            dataset["angle_to_next"].append(ftuv.vec_angle(dir1, dirN))
            dataset["distance_to_next"].append(ftuv.vec_distance(pos1, posN))
            dataset["next_stem"].append(next_stem)
        else:
            dataset["angle_to_next"].append("")
            dataset["distance_to_next"].append("")
            dataset["next_stem"].append("")
        dataset["outside_pk"].append(outside_pk)



def main():
    parser = fuc.get_rna_input_parser("Find pseudoknots in RNA structures, "
                                      "classify them into shapes and analyze "
                                      "their 3D architecturre.", "+",
                                      parser_kwargs={"conflict_handler":"resolve"})
    parser.add_argument("--pseudoknots", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--outfile-mode", choices=["w","a"], default='w', help="Overwrite ('w') or append ('a') to output file")
    parser.add_argument("--minlength", type= int, help= "Minimum length of each stem. "
                            "Stems with fewer base-pairs are treated as unpaired.",
                            default = 2)

    args =  parser.parse_args()
    args.pseudoknots=True
    rnas, filenames = fuc.cgs_from_args(args, rna_type="any",return_filenames = True, skip_errors=True)


    #variables for statistics
    unfold = 0
    without_pk = 0
    with_pk = 0
    pseudoknot_dataset = []
    pseudoknot_dataset_extended = defaultdict(list)

    for pos, rna in enumerate(rnas):
        try:
            without_pk, with_pk, unfold, pseudoknots = identification_pseudoknot\
                                            (rna, without_pk, with_pk, unfold,args.minlength)

            #count types of pseudoknots
            total_pk_g1 = len(pseudoknots)

            pk_sortclasses = {}
            other = 0
            other_pk =  []

            pk_classes, other = classify_pseudoknots(pseudoknots)

            print("pk_classes:")
            for key, pks in sorted(pk_classes.items()):
                print("{:<17s} {}\t{}".format(key, len(pks), ", ".join(map(str, pks))))
            print("other:"+"\t"+", ".join(map(str,other)))

            filename = str(filenames[pos]).split("/")[-1]

            entry={}
            for key, pks in pk_classes.items():
                entry[key] = len(pks)
            entry["other"] = len(other)
            entry["PK_other_structures"] = ",".join(map(str, other))
            entry["filename"] = filename
            entry["rnaname"] = rna.name
            print(filename)
            pseudoknot_dataset.append(entry)
            pk_id = 0
            for key, pks in pk_classes.items():
                for pk in pks:
                    pk_id+=1
                    extend_pk_description(pseudoknot_dataset_extended, filename,
                                          key, rna, pk, pk_id)


            for pk in other:
                pk_id+=1
                extend_pk_description(pseudoknot_dataset_extended, filename,
                                      "other", rna, pk, pk_id)
        except GraphIntegrityError:
            log.exception("Ignoring RNA %s; GraphIntegrityError", rna.name)
        except GraphConstructionError:
            log.exception("Ignoring RNA %s; GraphConstructionError", rna.name)
        except Exception:
            log.error("Error processing %s", rna.name)
            raise
    df1 = pandas.DataFrame(pseudoknot_dataset)
    df1.to_csv("pseudoknot_identification_genus2.csv", mode=args.outfile_mode, header=args.outfile_mode!="a", sep="\t")

    df2 = pandas.DataFrame(pseudoknot_dataset_extended)
    df2.to_csv("pseudoknot_identification_extended_genus2.csv", mode=args.outfile_mode, header=args.outfile_mode!="a", sep="\t")


    print("Structures with Pseudoknots: {}".format(with_pk))
    print("Structures without Pseudoknots: {}".format(without_pk))
    print("Structures unfold: {}".format(unfold))

if __name__ == "__main__":
    main()
