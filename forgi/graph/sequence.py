import logging
from collections import defaultdict
from . import residue as fgr

logger = logging.getLogger()

class _Smallest(object):
    """
    Smaller than everything, regardless of type
    """
    def __lt__(self, other):
        return True

class _Biggest(object):
    """
    Bigger than everything, regardless of type.
    """
    def __gt__(self, other):
        return True

def to_0_based(key):
    logger.debug("Converting key %s", key)
    if isinstance(key, slice):
        # Step
        step = key.step
        if step is None:
            step=1
        if abs(step)!=1:
            raise IndexError("Slices with step != 1 are not supported!")
        # Start
        start = key.start
        if start is None:
            start = None
        elif start > 0:
            start = start-1
        else:
            raise IndexError("Slices with start <1 are currently not defined.")
        # Stop
        if step<0:
            if key.stop is None:
                stop=None
            elif key.stop<-1:
                raise IndexError("Slices with negative stop and negative step are currently not defined.")
            else:
                stop = key.stop-2
        else:
            stop = key.stop
        key = slice(start, stop, step)
    elif isinstance(key, int):
        if key==0:
            raise IndexError("Index 0 not allowed in 1-based indexing")
        elif key>0:
            key=key-1
    else:
        raise TypeError("Invalid index type {}".format(type(key).__name__))
    logger.debug("Returning key %s", key)
    return key



def _resid_key(x):
    if x.resid[2] is None:
        return (x.resid[1], " ")
    else:
        return (x.resid[1], x.resid[2])

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
        reslist.sort(key=_resid_key)
    return chain_to_residues, resid_to_nucleotide

class _WMIndexer(object):
    def __init__(self, parent):
        self.parent = parent
    def __getitem__(self, key):
        return self.parent._getitem(key, True)

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
        self._seq = seq # A simple, 0-based string
        self._seqids = seqids
        mr, mnts = _sorted_missing_residues(missing_residues)
        self._missing_residues = mr
        self._missing_nts =  mnts

    def __getitem__(self, key):
        """
        Never return missing residues!
        """
        return self._getitem(key, False)

    def _getitem(self, key, include_missing):
        if isinstance(key, int):
            key=to_0_based(key)
            return self._seq[key]
        elif isinstance(key, fgr.RESID):
            try:
                i = self._seqids.index(key)
            except ValueError:
                if key in self._missing_nts:
                    nt = self._missing_nts[key]
                    if include_missing==False:
                        raise IndexError("No structure available for nucleotide '{}'."
                                         "For look-up including missing residues, use"
                                         "`.with_missing[key]`".format(key))
                    return nt
                raise IndexError("Nucleotide {} is not part of this RNA".format(key))
            else:
                return self._seq[i]
        elif isinstance(key, slice):
            return self._getslice(key, include_missing)
        else:
            raise TypeError("Wrong index type: {}".format(type(key).__name__))

    def _getslice(self, key, include_missing):
        """
        Do type checks and delegate to self._resid_slice or self._integer_slice
        """
        if key.step is not None and key.step not in [1, -1]:
            raise IndexError("The step in the slice has to be 1 or -1")
        if not isinstance(key.start, (type(None), int, fgr.RESID)):
            raise TypeError("Slice-start must be either an integer, a RESID or None")
        if not isinstance(key.stop, (type(None), int, fgr.RESID)):
            raise TypeError("Slice-stop must be either an integer, a RESID or None")
        if isinstance(key.start, fgr.RESID) or (key.start is None and isinstance(key.stop, fgr.RESID)):
            if not isinstance(key.stop, (fgr.RESID, type(None))):
                raise TypeError("Mixing integers and RESID for slicing is not allowed")
            return self._resid_slice(key, include_missing)
        return self._integer_slice(to_0_based(key), include_missing)

    def _integer_slice(self, key, include_missing):
        """
        This provides an implementation for include_missing=False
        and delegates to resid_slice otherwise.

        :param key: Already 0-based, no None in slice
        """
        if include_missing==False or not self._missing_nts:
            return self._seq[key]
        if key.start is None:
            startRES=None
        else:
            startRES = self._seqids[key.start]
        if key.stop is None:
            stopRES =  None
        else:
            stopRES  = self._seqids[key.stop]
        return self._resid_slice(slice(startRES, stopRES, key.step),  include_missing)

    def _resid_slice(self, key, include_missing):
        """
        This provides an implementation for inclue_missing=True
        and delegates to integer_slice otherwise.

        :param key: A validated resid slice.
        """
        left_res = None
        # Flag used later if missing_residues === True
        found_start=True
        if key.start is not None:
            try:
                start_i = self._seqids.index(key.start)
            except ValueError:
                if not include_missing:
                    raise IndexError("RESID {} is not present or no structure"
                                     "is available for it.".format(key.start))
                else:
                    start_i = None
                    if key.step!=-1:
                        found_start = False
        if key.start is None or start_i is None:
            if key.step==-1:
                start_i = len(self._seqids)
            else:
                start_i=0
        if key.stop is not None:
            try:
                stop_i = self._seqids.index(key.stop)
            except ValueError:
                if not include_missing:
                    raise IndexError("RESID {} is not present or no structure"
                                     "is available for it.".format(key.start))
                else:
                    stop_i = None
                    if key.step==-1:
                        found_start=False
        if key.stop is None or stop_i is None:
            if key.step==-1:
                stop_i=0
            else:
                stop_i = len(self._seqids)
        if not include_missing:
            # Convert to integers and return faster integer slice
            if key.step!=-1:
                stop_i+=1
            else:
                stop_i-=1
                if stop_i==-1:
                    stop_i=None
            key=slice(start_i, stop_i, key.step)
            logger.debug("Resid-slice mapped to integer-slice {}".format(key))
            return self._integer_slice(key, include_missing)
        else:
            seq = ""
            left_res = None
            start, stop = key.start, key.stop
            if key.step==-1:
                start_i, stop_i = stop_i, start_i
                start, stop = stop,start
            if start_i>0:
                left_res = self._seqids[start_i-1]
            logger.debug("Range is %s-%s", start_i, stop_i)
            for i in range(start_i, stop_i+1):
                try:
                    right_res = self._seqids[i]
                except IndexError:
                    right_res = None
                mr, mr_keys = self._missing_residues_between(left_res, right_res)
                for j,r in enumerate(mr_keys):
                    if not found_start and r==start:
                        logger.debug("Found start for i=%d, %s", i, r)
                        found_start=True
                    if found_start:
                        seq+=mr[j]
                    logger.debug("Comparing %s and %s", r, stop)
                    if r==stop:
                        logger.debug("Stop found! Seq=%s, step=%s", seq, key.step)
                        if key.step==-1:
                            seq="".join(reversed(seq))
                            logger.debug("Reversed to %s", seq)
                        return seq
                # If key.start and key.stop resp. are not modified residues,
                # start_i and stop_i are already accurrate, making
                # additional checks here unneccessary.
                if found_start and i<len(self._seq):
                    seq+=self._seq[i]
                left_res = right_res
                logger.debug("seq now is %s", seq)
            if key.step==-1:
                seq="".join(reversed(seq))
                logger.debug("Reversed to %s", seq)
            return seq

    def _missing_residues_between(self, from_resid, to_resid):
        if from_resid is not None and to_resid is not None:
            if from_resid.chain != to_resid.chain:
                left = self._missing_residues_between(from_resid, None)
                right = self._missing_residues_between(None, to_resid)
                return left[0]+right[0], left[1]+right[1]
        if from_resid is not None:
            from_key = _resid_key(from_resid)
            chain = from_resid.chain
        else:
            from_key=_Smallest()
        if to_resid is not None:
            to_key = _resid_key(to_resid)
            chain = to_resid.chain
        else:
            to_key=_Biggest()
        seq = ""
        res_ids = []
        for res1 in self._missing_residues[chain]:
            if _resid_key(res1)>from_key and _resid_key(res1)<to_key:
                seq+=self._missing_nts[res1]
                res_ids.append(res1)
        return seq, res_ids

    def to_resid(self, i):
        return self._seqids[i-1]

    def to_integer(self, resid):
        return self._seqids.index(resid)+1

    @property
    def with_missing(self):
        """
        Allows slices with missing residues!
        """
        if self._seqids:
            return _WMIndexer(self)
        else: # For performance reasons
            return self
