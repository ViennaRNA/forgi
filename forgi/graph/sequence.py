
def _sorted_missing_residues(list_of_dicts):
    chain_to_residues = defaultdict(list)
    resid_to_nucleotide = {}
    for res_dict in list_of_dicts:
        if res_dict["insertion"] is None:
            insertion=" "
        else:
            insertion = res_dict["insertion"]
        resid = fgr.RESID(res_dict["chain"], (' ', res_dict["ssseq"], insertion))
        chain_to_residues[res_dict["chain"]].append(resid)
        resid_to_nucleotide[resid] = res_dict["res_name"]
    for reslist in chain_to_residues.values():
        reslist.sort(key=lambda x: (x.resid[1], " ") if x.resid[2] is None
                                   else (x.resid[1], x.resid[2]))
    return chain_to_residues, resid_to_nucleotide

class Sequence(object):
    """
    There are two indexing conventions:
    *) 1-based indexing into the sequence for which structure information is available
       This uses indices of type integer.
       Integer based indexing can not address missing residues!
    *) Indexing using PDB-style numbers (They can start at any positive or negative number and
       may contain insertion codes or gaps.)
       This uses indices of type fgr.RESID
       PDB-Style indices may address missing residues (Reidues for which only sequence
       but no structure information is present).
    """
    def __init__(self, seq, seqids, missing_residues):
        """
        :param missing_residues: A list of dictionaries with the following keys:
                "res_name", "chain", "ssseq", "insertion"
        """
        self._seq = seq
        self._seqids = seqids
        self._missing_residues, self.missing_nts = _sorted_missing_residues(missing_residues)

    def __getitem__(self, key): 
        pass

    def to_resid(self, i):
        pass

    def to_integer(self, resid):
        pass

    def iter_seq(self, from_=None, to_=None, include_missing=False):
        pass
