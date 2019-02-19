"""
This module has functions which take a Bulge-Graph as input
and return a modified copy of it.

This is implemented as a class which is accessible as
BulgeGraph.transformed. This way we can properly inherit in CoarseGrainRNA
"""
import copy
import logging

from .sequence import MissingResidue, Sequence


log = logging.getLogger(__name__)


class _GCDummy(object):
    """
    Can be used in place of GraphConstruction for BG-Initialization
    """

    def __init__(self, defines, edges):
        self.defines = defines
        self.edges = edges


class BGTransformer(object):
    def __init__(self, bg):
        self.bg = bg

    def _without_elements(self, elems):
        """
        Return a copy of the BulgeGraph without the elements in elem.
        Their residues will be converted to missing residues.

        :param elems: A list of element names, e.g. ["s1","s2"]
        """
        raise NotImplementedError("This is still work in progress.")
        # We use the PDB numbering in new_defines, so we do not have to adjust indices if we remove
        # residues in the middle.
        resid_defines = {}
        for k, v in self.bg.defines:
            resid_defines[k] = list(map(self.bg.seq.to_resid, v))
        new_edges = copy.deepcopy(self.bg.edges)
        to_missing = []
        for elem in elems:
            if elem[0] != "i":
                raise NotImplementedError("TODO")
            else:
                stem1, stem2 = new_edges[elem]
                elem_define_a = list(
                    map(self.bg.seq.to_resid, self.bg.define_a(elem)))
                # Remove the iloop
                to_missing.extend(
                    self.bg.define_residue_num_iterator(elem, seq_ids=True))
                new_edges[stem1].remove(elem)
                new_edges[stem2].remove(elem)
                del new_edges[elem]
                # Merge stem2 into stem1
                new_edges[stem1] |= new_edges[stem2]
                del new_edges[stem2]
                if resid_defines[stem1][1] in elem_define_a:
                    # stem1 -IL - stem2
                    resid_defines[stem1][1] = resid_defines[stem2][1]
                    resid_defines[stem1][2] = resid_defines[stem2][2]
                else:
                    # stem2 - il - stem1
                    resid_defines[stem1][0] = resid_defines[stem2][0]
                    resid_defines[stem1][3] = resid_defines[stem2][3]
                del resid_defines[stem2]
                del resid_defines[elem]

    def condensed(self):
        """
        Return a condensed copy of the BulgeGraph.

        In the condensed BulgeGraph only the first (most 5-prime) nucleotide
        or base-pair of each element is retained, and the other nts/ base-pairs are
        converted to missing residues.
        In basepairs the first basepair contains the most 5' and most 3'
        nucleotide of the stem.
        """
        log.debug("Condensing BG with break-points %s",
                  self.bg.backbone_breaks_after)
        log.info("Condensing Graph %s", self.bg.to_dotbracket_string())
        new_defines = {}
        new_seqids = []
        new_seq = ""
        new_missing = self.bg.seq.with_missing.export_missing()
        new_i = 1
        for elem in self.bg.iter_elements_along_backbone():
            if not self.bg.defines[elem]:
                new_defines[elem] = []
            else:
                if elem in new_defines:  # Backwards strand for stems:
                    if len(new_defines[elem]) != 2:
                        log.error("%s", self.bg.edges[elem])
                    assert len(new_defines[elem]) == 2, "{} doesn't have len 2".format(
                        new_defines[elem])
                    assert elem[0] in "si"
                    fr, to = self.bg.defines[elem][2:]
                    if elem[0] == "s":
                        # to keep basepairing consistent, we keep the last (not the first)
                        # nt.
                        keep_i = to
                    else:
                        # To keep interiour loops consistent with the case where
                        # it is only at the second strand, keep the first nt here
                        keep_i = fr
                    new_defines[elem].extend([new_i, new_i])
                    log.debug("Extended new_defines for %s to %s",
                              elem, new_defines[elem])
                else:
                    fr, to = self.bg.defines[elem][:2]
                    keep_i = fr
                    new_defines[elem] = [new_i, new_i]
                    log.debug("Set new_defines for %s to %s",
                              elem, new_defines[elem])
                new_i += 1
                new_seq += self.bg.seq[keep_i]
                new_seqids.append(self.bg.seq.to_resid(keep_i))
                for i in range(fr, to + 1):
                    if i != keep_i:
                        seq_id = self.bg.seq.to_resid(i)
                        new_missing.append(
                            MissingResidue(seq_id, self.bg.seq[i]))
                    if i in self.bg.backbone_breaks_after:
                        if i >= keep_i:
                            new_seq += "&"
                        else:
                            new_seq = new_seq[:-1] + "&" + new_seq[-1]
                        log.debug("Breakpoint %s: new_seq now %s", i, new_seq)
        log.info("Condensing iteration done. Now creating condensed BG")
        graph_constr = _GCDummy(new_defines, copy.deepcopy(self.bg.edges))
        seq = Sequence(new_seq, new_seqids, new_missing,
                       self.bg.seq._modifications)
        return type(self.bg)(graph_constr, seq, name=self.bg.name + "_condensed", _dont_split=True)
