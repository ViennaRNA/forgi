from __future__ import unicode_literals


import logging
from collections import defaultdict
from . import residue as fgr

log = logging.getLogger()

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
    log.debug("Converting key %s", key)
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
    log.debug("Returning key %s", key)
    return key


def _insert_breakpoints_simple(seq, breakpoints, start=0, reverse = False):
    """
    :param seq: A sequence without any missing residues
    :param breakpoints: 0-based break points
    :param stop: The coordinate "on the forward strand",
                 i.e. the lower number of start and stop.
    """
    log.debug("Inserting breakpoints into %s", seq)
    breakpoints = sorted(breakpoints)
    if not breakpoints:
        return seq
    if start is None:
        start = 0
    out = []
    if reverse:
        seq=seq[::-1]
    log.debug("Inserting breakpoints into %s, start=%s, "
              "reversed=%s", seq, start, reverse)
    oldbp = 0
    for bp in breakpoints:
        log.debug("BP=%s", bp)
        bp-=start
        log.debug("bp-=start-->>>%s", bp)
        out.append(seq[oldbp:bp+1])
        oldbp=bp+1
    out.append(seq[bp+1:])
    out = [o for o in out if o]
    seq = "&".join(out)
    if reverse:
        seq=seq[::-1]
    log.debug("seq with break points: %s", seq)
    return seq

def _resid_key(x):
    if x.resid[2] is None:
        return (x.resid[1], " ")
    else:
        return (x.resid[1], x.resid[2])

def _sorted_missing_residues(list_of_dicts):
    chain_to_residues = defaultdict(list)
    resid_to_nucleotide = {}
    for res_dict in list_of_dicts:
        if res_dict["model"] not in [None, 1, "A"]:
            continue
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
    def __len__(self):
        return len(self.parent)+len(self.parent._missing_nts)
    def update_dotbracket(self, struct):
        """
        Given a dotbracket_string of the same length as the sequence,
        update it with '-' for missing residues and '&' for breakpoints.

        :returns: A string
        """
        return self.parent._missing_into_db(struct)

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
    def __init__(self, seq, seqids, missing_residues=[]):
        """
        :param missing_residues: A list of dictionaries with the following keys:
                "res_name", "chain", "ssseq", "insertion"
        """
        log.debug("Sequence initialized with %s, %s, %s", seq, list(map(fgr.resid_to_str, seqids)), missing_residues)
        # Uses 0-based indexing
        i=0
        self._breaks_after = []
        remaining = seq
        while remaining:
            c, remaining = remaining[0], remaining[1:]
            if c=='&':
                self._breaks_after.append(i-1)
            else:
                i+=1
        log.debug("Break-points for seq %s are: %s", seq, self._breaks_after)
        self._seq = seq.replace('&', '')
        self._seqids = seqids
        mr, mnts = _sorted_missing_residues(missing_residues)
        self._missing_residues = mr
        self._missing_nts =  mnts # A dict seq_id:nt

    def __str__(self):
        return self[:]

    def __repr__(self):
        try:
            qname = type(self).__qualname__
        except AttributeError:
            qname = type(self).__name__
        return "<{} object with ._seq={}>".format(qname, str(self))

    def __bool__(self):
        return bool(self._seq)
    __nonzero__ = __bool__

    @property
    def backbone_breaks_after(self):
        """
        Breakpoints as 1-based indices
        """
        return [ bp+1 for bp in self._breaks_after ]

    def is_valid(self):
        wrong_chars = set(self.seq)-set("AUGCaugc")
        if wrong_chars:
            log.info("Illegal characters are {}".format(wrong_chars))
            return False
        return True

    def __getitem__(self, key):
        """
        Never return missing residues!
        """
        return self._getitem(key, False)

    def __len__(self):
        return len(self._seq)

    def __eq__(self, other):
        log.debug("{} =?= {}".format(repr(self), repr(other)))
        log.debug("type of other=%s, type of string literal in this module is %s", type(other).__name__, type("").__name__)
        if isinstance(other, type(self)):
            return (self._seq == other._seq and
                    self._seqids == other._seqids and
                    self._missing_nts == other._missing_nts and
                    self._breaks_after == other._breaks_after)
        elif isinstance(other, type("")):
            return str(self)==other
        else:
            return NotImplemented

    def _getitem(self, key, include_missing):
        log.debug("_getitem called for %s, include_missing=%s", key, include_missing)
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
            seq = self._seq[key]
            if key.step==-1:
                start=key.stop
                if start is not None:
                    start+=1
            else:
                start=key.start
            return _insert_breakpoints_simple(seq, self._breaks_after, start, key.step==-1)
        if key.start is None:
            startRES=None
        else:
            startRES = self._seqids[key.start]
        if key.stop is None:
            stopRES =  None
        else:
            if key.step==-1:
                stopRES  = self._seqids[key.stop+1]
            else:
                if key.stop>0:
                    stopRES  = self._seqids[key.stop-1]
                else:
                    stopRES  = self._seqids[key.stop]
        key2 = slice(startRES, stopRES, key.step)
        log.debug("Converted integer_slice %s to resid-slice %s", key, key2)
        return self._resid_slice_with_missing(key2)

    def _resid_slice_without_missing(self, key):
        if key.step==-1:
            default_start = len(self._seqids)
            default_stop = None
        else:
            default_start = 0
            default_stop = len(self._seqids)
        start_i, _ = self._resid_to_index(key.start, default_start, False)
        stop_i, _ = self._resid_to_index(key.stop, default_stop, False)
        if key.step==-1:
            if stop_i is not None:
                stop_i-=1
        else:
            stop_i+=1
        key=slice(start_i, stop_i, key.step)
        log.debug("Resid-slice mapped to integer-slice {}".format(key))
        return self._integer_slice(key, False)


    def _resid_to_index(self, resid, default, default_on_error=False):
        """
        :returns: a tuple index, error_raised
        """
        if resid is None:
            return default, False
        try:
            return self._seqids.index(resid), False
        except ValueError:
            if default_on_error:
                return default, True
            else:
                raise IndexError("RESID {} is not present or no structure"
                                 "is available for it.".format(resid))


    def _resid_slice(self, key, include_missing):
        """
        This provides an implementation for inclue_missing=True
        and delegates to integer_slice otherwise.

        :param key: A validated resid slice.
        """
        if include_missing:
            return self._resid_slice_with_missing(key)
        else:
            return self._resid_slice_without_missing(key)

    def _resid_slice_with_missing(self, key):
        if key.step==-1:
            start, stop = key.stop, key.start
        else:
            start, stop = key.start, key.stop
        seqids = list(self._iter_resids_with_missing(start, stop))
        if key.step==-1:
            seqids.reverse()
        seq = ""
        for seqid in seqids:
            if seqid == "&":
                seq+="&"
            else:
                try:
                    seq+=self._seq[self._seqids.index(seqid)]
                except ValueError:
                    seq+=self._missing_nts[seqid]
        return "".join(seq)

    def _missing_into_db(self, dotbracket):
        db = ""
        check_i = 0
        for seqid in self._iter_resids_with_missing(None, None):
            if seqid == "&":
                db+="&"
            else:
                try:
                    i = self._seqids.index(seqid)
                    assert i==check_i
                except ValueError:
                    db+="-"
                else:
                    db+=dotbracket[i]
                    check_i+=1
        return "".join(db)


    def _iter_resids_with_missing(self, start, stop):
        left_res = None
        start_i, start_is_missing = self._resid_to_index(start, 0, True)
        stop_i, stop_is_missing   = self._resid_to_index(stop, len(self._seqids), True)
        if not start_is_missing and start is not None:
            left_res = start

        # Flag used later if missing_residues == True
        found_start=not start_is_missing
        seq = ""
        if left_res is not None:
            yield left_res
            #seq+=self._seq[self._seqids.index(left_res)]
            #log.debug("Added first residue %s: %s", left_res, seq)
            start_i+=1
        log.debug("Iterating over sequence %s-%s", start_i, stop_i+1)
        for i in range(start_i, stop_i+1):
            try:
                right_res = self._seqids[i]
            except IndexError:
                right_res = None
            # Some consistency checks:
            if left_res is not None and right_res is not None and left_res.chain!=right_res.chain:
                if i-1 not in self._breaks_after:
                    raise ValueError("Inconsist break-points discovered."
                                     "Missing breakpoint at position {} "
                                     "between chain {} and chain {}".format(i-1,
                                     left_res.chain,right_res.chain))
            if i-1 in self._breaks_after:
                if left_res is not None and right_res is not None and left_res.chain==right_res.chain:
                    raise ValueError("Inconsist break-points discovered."
                                     "Breakpoint at position {} "
                                     "in the middle of chain {}".format(i-1,
                                      left_res.chain))
            mr, mr_keys = self._missing_residues_between(left_res, right_res)
            old_r = left_res
            for j,r in enumerate(mr_keys):
                if not found_start and r==start:
                    log.debug("Found start for i=%d, %s", i, r)
                    found_start=True
                    old_r = None
                log.debug("Now processing missing residue %s", r)
                if found_start:
                    if old_r is not None and old_r.chain != r.chain:
                        yield "&"
                        log.debug("Breakpoint inserted!")
                    log.debug("Adding missing residue %s", mr[j])
                    yield r
                log.debug("Comparing %s and stop=%s", r, stop)
                if r==stop:
                    log.debug("Stop found! Seq=%s", seq)
                    return
                old_r = r
            # If key.start and key.stop resp. are not modified residues,
            # start_i and stop_i are already accurrate, making
            # additional checks here unneccessary.
            if found_start and i<len(self._seq):
                log.debug("Adding regular residue %s", self._seq[i])
                if mr_keys: #No missing residues
                    if r.chain!=right_res.chain:
                        yield "&"
                        log.debug("Breakpoint inserted!")
                else:
                    if left_res is not None and right_res is not None and left_res.chain!=right_res.chain:
                        yield "&"
                        log.debug("Breakpoint inserted!")
                yield right_res
            left_res = right_res
            log.debug("seq now is %s", seq)
        # If start or stop are not part of sequence,
        # make sure they were seen in missing residues
        if not found_start or stop_is_missing:
            raise IndexError("At least one of the residues {} and {} "
                             "is not part of this RNA".format(key.start, key.stop))
        return

    def _missing_residues_between(self, from_resid, to_resid):
        log.debug("Searching missing residues between %s  and %s", from_resid, to_resid)
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
