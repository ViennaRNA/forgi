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

    def without_elements( self, elems ):
        """
        Return a copy of the BulgeGraph without the elements in elem.
        Thier residues will be converted to missing residues.

        :param elems: A list of element names, e.g. ["s1","s2"]
        """
        raise NotImplementedError("TODO")

    def condensed(self):
        """
        Return a copy of the bulge Graph, where for each CoarseGrained Element,
        only the first nucleotide or base-pair is retained, and the other
        nts/ base-pairs are converted to missing residues.
        """
        log.debug("Condensing BG with break-points %s", self.bg.backbone_breaks_after)
        new_defines = {}
        new_seqids = []
        new_seq = ""
        new_missing = self.bg.seq.with_missing.export_missing()
        new_i = 1
        for elem in self.bg.iter_elements_along_backbone():
            if not self.bg.defines[elem]:
                new_defines[elem]=[]
            else:
                if elem in new_defines: # Backwards strand for stems:
                    assert len(new_defines[elem])==2
                    assert elem[0] in "si"
                    fr, to = self.bg.defines[elem][2:]
                    if elem[0]=="s":
                        # to keep basepairing consistent, we keep the last (not the first)
                        # nt.
                        keep_i = to
                    else:
                        # To keep interiour loops consistent with the case where
                        # it is only at the second strand, keep the first nt here
                        keep_i = fr
                    new_defines[elem].extend([new_i, new_i])
                else:
                    fr, to = self.bg.defines[elem][:2]
                    keep_i = fr
                    new_defines[elem]=[new_i, new_i]
                new_i+=1
                new_seq+=self.bg.seq[keep_i]
                new_seqids.append(self.bg.seq.to_resid(keep_i))
                for i in range(fr, to+1):
                    if i != keep_i:
                        seq_id = self.bg.seq.to_resid(i)
                        new_missing.append(MissingResidue(seq_id, self.bg.seq[i]))
                    if i in self.bg.backbone_breaks_after:
                        if i>keep_i:
                            new_seq+="&"
                        else:
                            new_seq=new_seq[:-1]+"&"+new_seq[-1]
        log.info("Condensing iteration done. Now creating condensed BG")
        graph_constr = _GCDummy(new_defines, copy.deepcopy(self.bg.edges))
        seq = Sequence(new_seq, new_seqids, new_missing, self.bg.seq._modifications)
        return type(self.bg)(graph_constr, seq, name=self.bg.name+"_condensed", _dont_split=True)
