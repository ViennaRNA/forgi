from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import object
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from pprint import pprint

import json
import warnings
import sys
from collections import Counter, defaultdict, namedtuple
import itertools as it
import logging

from logging_exceptions import log_to_exception

from forgi.graph.bulge_graph import RESID

log = logging.getLogger(__name__)


class DSSRLookupError(LookupError):
    pass


class WrongChain(LookupError):
    pass


def dssr_to_pdb_resid(dssr_resid):
    if ":" in dssr_resid:
        model, dssr_resid = dssr_resid.split(":")
        if model != "1":
            return None
    if '.' in dssr_resid:
        chain, _, resid = dssr_resid.partition(".")
    else:
        chain = None
        resid = dssr_resid
    if "^" in resid:
        resid, _, letter = resid.partition("^")
    else:
        letter = " "
    nat = []
    for c in reversed(resid):
        if c.isdigit() or c == "-":
            nat.append(c)
        else:
            break
    resid = "".join(reversed(nat))
    return RESID(chain, (" ", int(resid), letter))


class DSSRAnnotation(object):
    def __init__(self, dssr, cg):
        """
        Store annotation from dssr associated with a coarse grained RNA.

        For DSSR, see: doi: 10.1093/nar/gkv716

        .. note::
            This class is currently compatible with v1.5.7-2016jun16 of dssr.

        :param dssr: A *.json file returned by `x3dna-dssr --json`, or a string holding the json
                     or a dictionary containing the parsed json.
        :param cg: The CoarseGrainedRNA that corresponds to
                   the same pdb file x3dna-dssr was run on.
        """
        if isinstance(dssr, type("")):
            if dssr[0] != "{":
                with open(dssr) as f:
                    dssr = f.read()
            dssr = json.loads(dssr)
        elif not isinstance(dssr, dict):
            raise TypeError(
                "dssr must be an string- or dict instance, not {}".format(type(dssr)))
        self._dssr = dssr
        self._cg = cg

    def aminor_interactions(self):
        if "Aminors" not in self._dssr:
            return
        for aminor in self._dssr["Aminors"]:
            nt1,_,nt23 = aminor["desc_long"].split()
            nt2, nt3 = map(dssr_to_pdb_resid, nt23.split(","))
            nt1 = dssr_to_pdb_resid(nt1)
            if nt1 not in self._cg.seq or nt2 not in self._cg.seq or nt3 not in self._cg.seq:
                log.debug("Ignoring A-Minor interaction outside this connected component: %s->%s-%s", nt1, nt2, nt3)
                continue
            stem = self._cg.get_elem(nt2)
            if stem!=self._cg.get_elem(nt3):
                log.warning("DSSR assumes {} and {} are in a stem, but we assume they are not. Ignoring this AMinor annotation".format(nt2, nt3))
                continue
            yield nt1, stem


    def noncanonical_pairs(self, elem=None):
        """
        Yield all nucleotides which are in a non-canonical basepair
        or a canonical basepair outside a stem.

        This ignores interactions with other RNA chains that are not
        part of the CoarseGrainRNA. Thus, yielded nucleotides are valid indices
        to the CoarseGrainRNA's `.seq` object.

        :param elem: An element name, e.g. "i0".
                     If given, only find basepairs within this secondary structure element.
        :param yields: A triple: `nt1, nt2, basepair-type`
                       nt1 and nt2 are forgi.graph.residue.RESID instances (PDB numbering scheme).
                       The basepair-type is given in the LW notation,
                       e.g. cWH for cis basepairing between the Watson-Crick anfd the Hoogsteen edge.
                       or   tHS for trans basepairing between the Hoogsteen and the Sugar edge.
        """
        for pair in self._dssr["pairs"]:
            nt1 = dssr_to_pdb_resid(pair["nt1"])
            nt2 = dssr_to_pdb_resid(pair["nt2"])
            if nt1 is None or nt2 is None:
                continue
            bp_type = pair["LW"]
            try:
                s1 = self._cg.seq.with_missing[nt1]
                s2 = self._cg.seq.with_missing[nt2]
                log.debug("Basepair is %s -%s- %s", s1, bp_type, s2)
            except IndexError:
                continue
            if elem is None:
                if bp_type == "cWW":
                    partner = self._cg.pairing_partner(nt2)
                    if nt1 == partner:
                        continue
                    elif self._cg.get_node_from_residue_num(nt1)[0] == "s" and self._cg.get_node_from_residue_num(nt2)[0] == "s":
                        log.warning("Pairing partner of %s is %s, but "
                                    "DSSR has %s in interaction. Base-triple?", nt2, partner, nt1)
                yield nt1, nt2, bp_type
            else:
                try:
                    if self._cg.get_node_from_residue_num(nt1) == elem and self._cg.get_node_from_residue_num(nt2) == elem:
                        yield nt1, nt2, bp_type
                except ValueError:
                    pass

    def stacking_nts(self):
        """
        What nucleotides stack according to DSSR

        :returns: A list of tuples of stacking nts.
        """
        if "stacks" not in self._dssr:
            return []
        out = []
        for stack in self._dssr["stacks"]:
            nts = stack["nts_long"].split(",")
            resids = []
            for nt in nts:
                resid = dssr_to_pdb_resid(nt)
                if resid is None:
                    continue
                resids.append(resid)
            for i, resid in enumerate(resids):
                try:
                    out.append((resid, resids[i+1]))
                except IndexError:
                    pass
        return out

    def coaxial_stacks(self):
        """
        Get the element name in the CoarseGrainRNA for all coaxial stacks found by DSSR.
        """
        if "coaxStacks" not in self._dssr:
            return []
        cg_stacks = []
        for dssr_stack in self._dssr["coaxStacks"]:
            cg_stack = []
            try:
                for dssr_stem in dssr_stack["stem_indices"]:
                    cg_stack.append(self.cg_stem(dssr_stem))
            except WrongChain:
                continue
            except RuntimeError:
                continue
            else:
                log.debug("Coaxial stack %s mapped to forgi %s", dssr_stack, cg_stack)
                cg_stacks.append(cg_stack)
        log.debug("coaxial stacks are %s", cg_stacks)
        return cg_stacks

    def cg_stem(self, dssr_stem):
        """
        Get the stem define in the CoarseGrainRNA that corresponds to the stem id in the dssr-format.

        :param dssr_stem: INT the stem in the DSSR Anntotation.
        """
        log.debug("Mapping DSSR stem %s to forgi", dssr_stem)
        if "stems" not in self._dssr:
            raise DSSRLookupError("The DSSR object does not contain any stem!")
        for stem_obj in self._dssr["stems"]:
            if stem_obj["index"] == dssr_stem:
                break
        else:
            raise DSSRLookupError("No stem with index {}".format(dssr_stem))
        log.debug("Found stem %s&%s", stem_obj["strand1"], stem_obj["strand2"])

        cg_stems = Counter()  # See, if the dssr_stems maps to more than 1 cg-stem
        for pair in stem_obj["pairs"]:
            res1 = dssr_to_pdb_resid(pair["nt1"])
            res2 = dssr_to_pdb_resid(pair["nt2"])
            if res1 is None or res2 is None:
                continue
            log.debug("Contains pair %s-%s", res1, res2)
            if self._cg.chains and (res1.chain not in self._cg.chains or res2.chain not in self._cg.chains):
                e = WrongChain()
                with log_to_exception(log, e):
                    log.error("Wrong chain: res1={}, res2={}, cg.chains={}".format(
                        res1, res2, self._cg.chains))
                raise e
            i1 = self._cg.seq.to_integer(res1)
            i2 = self._cg.seq.to_integer(res2)
            nodes = self._cg.nucleotides_to_elements([i1, i2])
            for node in nodes:
                cg_stems[node] += 1
        if not cg_stems:
            raise RuntimeError(
                "No stem matching dssr_stem {}.".format(dssr_stem))
        most_common = cg_stems.most_common()
        if len(most_common) > 1:
            extra_info = ""
            for d in cg_stems.keys():
                if d[0] == "i":
                    extra_info += "\n{} is {}:".format(d,
                                                       self._cg.get_define_seq_str(d))
                    extra_info += "\n\t" + self._cg.seq + "\n\t" + \
                        self._cg.to_dotbracket_string() + "\n\t"
                    resnums = list(self._cg.define_residue_num_iterator(d))
                    for i in range(len(self._cg.seq)):
                        pos = i + 1
                        if pos in resnums:
                            extra_info += "^"
                        else:
                            extra_info += " "

            warnings.warn("dssr_stem {} maps to more than one cg element: {} {}".format(
                dssr_stem, list(cg_stems.keys()), extra_info))
        for mc in most_common:
            if mc[0][0] == "s":
                return mc[0]
        raise RuntimeError("No stem matching dssr_stem {}, only single stranded region: {}.".format(
            dssr_stem, list(cg_stems.keys())))

    def compare_coaxial_stack_annotation(self, forgi_method="Tyagi", allow_single_bp=False):
        """
        Compare the coaxial stack annotation between the DSSR method and the forgi method.
        :param forgi_method: "Tyagi" or "CG". Method for stack detection in forgi.
        """
        stacks_dssr = set()
        stacks_forgi = set()
        Stack = namedtuple('Stack', ['stems', 'forgi', 'dssr'])
        dssr_stacks = self.coaxial_stacks()
        for stack in dssr_stacks:
            if stack[0] == stack[1]:
                stacks_dssr.add(Stack(tuple(stack), "one helix", "stacking"))
            else:
                bulges = (self._cg.edges[stack[0]] & self._cg.edges[stack[1]])
                if not bulges:
                    stacks_dssr.add(
                        Stack(tuple(stack), "not connected", "stacking"))
                else:
                    bulge = bulges.pop()  # We only look at one bulge. If 2 stems are connected by more than 1 bulge, cg.is_stacking should give the same result in both cases. TODO: Write a test for tihis case!
                    if self._cg.is_stacking(bulge, forgi_method):
                        curr_stack = Stack(
                            tuple(stack), "stacking", "stacking")
                        stacks_dssr.add(curr_stack)
                        stacks_forgi.add(curr_stack)
                    else:
                        stacks_dssr.add(
                            Stack(tuple(stack), "not stacking", "stacking"))
        for d in self._cg.defines:
            if d[0] not in "mi":
                continue
            if self._cg.is_stacking(d, forgi_method):
                s1, s2 = self._cg.connections(d)
                if not allow_single_bp and (self._cg.stem_length(s1) == 1 or self._cg.stem_length(s2) == 1):
                    continue
                if [s1, s2] not in dssr_stacks and [s2, s1] not in dssr_stacks:
                    stacks_forgi.add(
                        Stack((s1, s2), "stacking", "not stacking"))
        return stacks_forgi, stacks_dssr

    def stacking_loops(self):
        """
        Returns a list of loop (multiloop and interior-loop) elements,
        along which the adjacent stems stack.

        Note that not all nucleotides of these elements have to participate in
        the stack and additional nucleotides can be used for the stack.

        I fact, this only checks whether one nucleotide of each of the first
        (seen from the side of the loop element) basepairs of the adjacent
        stems are in a continuouse stack.

        :returns: A list of element names
        """
        stacking = []
        for loop in self._cg.defines:
            if loop[0] not in "im":
                continue
            nts_elem = []
            define_a = self._cg.define_a(loop)
            pos1 = define_a[0]
            helix1 = [self._cg.seq.to_resid(pos1), self._cg.seq.to_resid(
                self._cg.pairing_partner(pos1))]
            pos2 = define_a[1]
            helix2 = [self._cg.seq.to_resid(pos2), self._cg.seq.to_resid(
                self._cg.pairing_partner(pos2))]
            for stack in self._dssr.get("stacks", []):
                nts_dssrstack = []
                for nt in stack["nts_long"].split(","):
                    nt = dssr_to_pdb_resid(nt)
                    if nt is None:
                        continue
                    nts_dssrstack.append(nt)
                if any(nt in nts_dssrstack for nt in helix1) and any(nt in nts_dssrstack for nt in helix2):
                    log.debug("elem %s is stacking ", loop)
                    stacking.append(loop)
                    break
        log.debug(stacking)
        return stacking

    def basepair_stacking(self, forgi_method="Tyagi"):
        # if "_" in self._cg.name:
        #    chainname="chain_"+self._cg.name.split("_")[-1]
        #    forgi_chain=self._cg.name.split("_")[-1]
        # else:
        if True:
            chainname = "all_chains"
            forgi_chain = None
        dssr_helices = []
        for stem in self._dssr.get("stems", []):
            dssr_helices.append(set())
            for pair in stem["pairs"]:
                nt1 = dssr_to_pdb_resid(pair["nt1"])
                nt2 = dssr_to_pdb_resid(pair["nt2"])
                if nt1 is None or nt2 is None:
                    continue
                if forgi_chain is None or nt1[0] == forgi_chain:
                    dssr_helices[-1].add(nt1)
                if forgi_chain is None or nt2[0] == forgi_chain:
                    dssr_helices[-1].add(nt2)
            if not dssr_helices[-1]:
                del dssr_helices[-1]
        for stack in self._dssr.get("stacks", []):
            dssr_helices.append(set())
            for nt in stack["nts_long"].split(","):
                nt1 = dssr_to_pdb_resid(nt)
                if nt1 is not None and forgi_chain is None or nt1[0] == forgi_chain:
                    dssr_helices[-1].add(nt1)
            if not dssr_helices[-1]:
                del dssr_helices[-1]
        # Merge stacks that overlap.
        while True:
            for i, j in it.combinations(range(len(dssr_helices)), 2):
                stack_bag1 = dssr_helices[i]
                stack_bag2 = dssr_helices[j]
                if stack_bag1 & stack_bag2:
                    stack_bag1 |= stack_bag2
                    del dssr_helices[j]
                    break
            else:
                break
        element_helices = self._cg.get_stacking_helices(forgi_method)
        forgi_helices = []
        for helix in element_helices:
            forgi_helices.append(set())
            for d in helix:
                forgi_helices[-1] |= set(
                    n for n in self._cg.define_residue_num_iterator(d, seq_ids=True))
        print("FH", forgi_helices)
        for i, d_helix in enumerate(dssr_helices):
            no_match = True
            for j, f_helix in enumerate(forgi_helices):
                if d_helix & f_helix:
                    no_match = False
                    print ("dssr helix {} matches forgi helix {}. Only forgi: {}nts, only dssr: {}nts, both: {}nts".format(
                        i, j, len(f_helix - d_helix), len(d_helix - f_helix), len(d_helix & f_helix)))
            if no_match:
                print("dssr helix {}: No match ({}nts)".format(i, len(d_helix)))
        for j, f_helix in enumerate(forgi_helices):
            no_match = True
            for i, d_helix in enumerate(dssr_helices):
                if d_helix & f_helix:
                    no_match = False
                    print ("forgi helix {} matches dssr helix {}. Only forgi: {}nts, only dssr: {}nts, both: {}nts".format(
                        j, i, len(f_helix - d_helix), len(d_helix - f_helix), len(d_helix & f_helix)))
            if no_match:
                print("forgi helix {}: No match ({}nts)".format(j, len(f_helix)))
        print("FH", forgi_helices)
        print ("seq   " + self._cg.seq)
        print ("forgi " + self._cg.to_dotbracket_string())
        display_f = []
        for helix in forgi_helices:
            display_f.append([self._cg.seq_id_to_pos(nt) - 1 for nt in helix])
        display_f.sort(key=lambda x: min(x))
        helixstri = "forgi "
        for i in range(len(self._cg.seq)):
            for j, dp in enumerate(display_f):
                if i in dp:
                    # Start with 0 (ord("0")==48 and display only characters between chr(33)="!" and chr(126)="~"
                    chrid = j + 15
                    if chrid > 93:
                        chrid = chrid % (93)
                    chrid = chrid + 33

                    helixstri += chr(chrid)
                    break
            else:
                helixstri += " "
        print(helixstri)
        print ("DSSR  " + self._dssr["dbn"][chainname]["sstr"])
        display_d = []
        # print(dssr_helices)
        for helix in dssr_helices:
            try:
                display_d.append(
                    [self._cg.seq_id_to_pos(nt) - 1 for nt in helix])
            except ValueError as e:
                log.exception("Exception occurred")
        display_d.sort(key=lambda x: min(x))
        helixstri = "dssr  "
        for i in range(len(self._cg.seq)):
            for j, dp in enumerate(display_d):
                if i in dp:
                    chrid = j + 15
                    if chrid > 93:
                        chrid = chrid % (93)
                    chrid = chrid + 33

                    helixstri += chr(chrid)
                    break
            else:
                helixstri += " "
        print(helixstri)
        print(dssr_helices)
        dssr_to_elems = [self._cg.nucleotides_to_elements(
            list(map(self._cg.seq.to_integer, helix))) for helix in dssr_helices]
        print(dssr_to_elems)

    def compare_dotbracket(self):
        if "_" in self._cg.name:
            chainname = "chain_" + self._cg.name.split("_")[-1]
        else:
            chainname = "all_chains"
        print ("seq   " + self._cg.seq)
        print ("forgi " + self._cg.to_dotbracket_string())
        print ("DSSR  " + self._dssr["dbn"][chainname]["sstr"])
