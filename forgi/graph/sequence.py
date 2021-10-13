from __future__ import absolute_import, unicode_literals, division
from __future__ import print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input, # pylint: disable=W0611
                      map, next, oct, pow, range, round, # pylint: disable=W0611
                      str, super, zip) # pylint: disable=W0611

try:
    from collections.abc import Sequence as SequenceABC
except ImportError:
    from collections import Sequence as SequenceABC

import logging
from collections import defaultdict, Counter
from string import ascii_lowercase, ascii_uppercase
from functools import partial
import inspect

from logging_exceptions import log_to_exception
from . import residue as fgr
log = logging.getLogger(__name__)

VALID_CHAINIDS = ascii_uppercase + ascii_lowercase


class MissingResidue(object):
    def __init__(self, resid, res_name):
        """
        :param resid: A fgr.RESID instance of a string.
        :param resname: A string. Probably one of "AUGC".
        """
        if not isinstance(resid, fgr.RESID):
            log.debug("Type of %s is %s", resid, type(resid).__name__)
            resid = fgr.resid_from_str(resid)
        self.resid = resid
        self.res_name = res_name

    def to_bg_string(self):
        return "missing {} {}".format(fgr.resid_to_str(self.resid), self.res_name)

    @classmethod
    def from_bg_fields(cls, parts):
        """
        Used during loading of the bg_file.

        :param parts: A list of strings. If the first is not "missing",
                      None is returned
        :returns: Either a MissingResidue instance or None
        """
        if parts[0] != "missing":
            return None
        return cls(fgr.resid_from_str(parts[1]), parts[2])


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
            step = 1
        if abs(step) != 1:
            raise IndexError("Slices with step != 1 are not supported!")
        # Start
        start = key.start
        if start is None:
            start = None
        elif start > 0:
            start = start - 1
        else:
            raise IndexError("Slices with start <1 are currently not defined.")
        # Stop
        if step < 0:
            if key.stop is None:
                stop = None
            elif key.stop < -1:
                raise IndexError(
                    "Slices with negative stop and negative step are currently not defined.")
            else:
                stop = key.stop - 2
        else:
            stop = key.stop
        key = slice(start, stop, step)
    elif isinstance(key, int):
        if key == 0:
            raise IndexError("Index 0 not allowed in 1-based indexing")
        elif key > 0:
            key = key - 1
    else:
        raise TypeError("Invalid index type {}".format(type(key).__name__))
    log.debug("Returning key %s", key)
    return key


def _insert_breakpoints_simple(seq, breakpoints, start=0, reverse=False):
    """
    :param seq: A sequence without any missing residues
    :param breakpoints: 0-based break points
    :param start: The coordinate "on the forward strand",
                 i.e. the lower number of start and stop.
    """
    log.debug("Inserting breakpoints %s into %s", seq, breakpoints)
    breakpoints = sorted(breakpoints)
    if not breakpoints:
        return seq
    if start is None:
        start = 0
    out = []
    if reverse:
        seq = seq[::-1]
    log.debug("Inserting breakpoints into %s, start=%s, "
              "reversed=%s", seq, start, reverse)
    oldbp = 0
    for bp in breakpoints:
        log.debug("BP=%s", bp)
        bp -= start
        log.debug("bp-=start-->>>%s", bp)
        out.append(seq[oldbp:bp + 1])
        oldbp = bp + 1
    out.append(seq[bp + 1:]) # pylint: disable=undefined-loop-variable
    out = [o for o in out if o]
    seq = "&".join(out)
    if reverse:
        seq = seq[::-1]
    log.debug("seq with break points: %s", seq)
    return seq


def _resid_key(x):
    if x.resid[2] is None:
        return (x.resid[1], " ")
    else:
        # Convert to newstring on python2 See: https://github.com/PythonCharmers/python-future/issues/110
        return (x.resid[1], str(x.resid[2]))


def _sorted_missing_residues(list_of_dicts):
    """
    :param list_of_dicts: A list of dicts, as in my PR to biopython
                          or with the two keys "RESID" and "res_name",
    """
    chain_to_residues = defaultdict(list)
    resid_to_nucleotide = {}
    for res_dict in list_of_dicts:
        log.debug("Processing %s", res_dict)
        if isinstance(res_dict, MissingResidue):
            resid = res_dict.resid
            chain = res_dict.resid.chain
            res_name = res_dict.res_name
        elif "RESID" in res_dict:
            resid = res_dict["RESID"]
            if not isinstance(resid, fgr.RESID):
                resid = fgr.resid_from_str(resid)
            chain = resid.chain
            res_name = res_dict["res_name"]
        else:
            if res_dict["model"] not in [None, 1, "A"]:
                log.info("Invalid missing residue %s", res_dict)
                continue
            if res_dict["insertion"] is None:
                insertion = " "
            else:
                insertion = res_dict["insertion"]
            chain = res_dict["chain"]
            resid = fgr.RESID(chain, (' ', res_dict["ssseq"], insertion))
            res_name = res_dict["res_name"]
        chain_to_residues[chain].append(resid)
        resid_to_nucleotide[resid] = res_name
    for reslist in chain_to_residues.values():
        log.debug("Sorting %s", reslist)
        reslist.sort(key=_resid_key)
    return chain_to_residues, resid_to_nucleotide


class _IndexHelper(object):
    # pylint: disable=protected-access
    @property
    def flag(self):
        raise NotImplementedError("Must be set in subclass")

    def __init__(self, parent):
        self.parent = parent

    def __getitem__(self, key):
        log.debug("%s.__getitem__ called. %s=True",
                  type(self).__name__, self.flag)
        return self.parent._getitem(key, **{self.flag: True})

    def _getitem(self, key, include_missing=False, show_modifications=False):
        kwargs = {"include_missing": include_missing,
                  "show_modifications": show_modifications}
        kwargs.update({self.flag: True})
        log.debug("%s._getitem called. flags: %s", type(self).__name__, kwargs)
        return self.parent._getitem(key, **kwargs)

    def __getattr__(self, attr):
        """
        If attr is a method that supports **kwargs or has an attribute with
        the same name as the value of self.flag, then we set this
        attributre to True
        """
        log.debug("Getattr called for %s", attr)
        f = getattr(self.parent, attr)
        if callable(f):
            argspec = inspect.getargspec(f)
            log.debug("For function %s: args are %s", f, argspec.args)
            if argspec.keywords or self.flag in argspec.args:
                log.debug("setting flag %s", self.flag)
                kwargs = {self.flag: True}
                f = partial(f, **kwargs)
        return f
    def __str__(self):
        return str(self[:])

class _WMIndexer(_IndexHelper):
    # pylint: disable=protected-access

    flag = "include_missing"

    def __len__(self):
        return len(self.parent) + len(self.parent._missing_nts)

    def update_dotbracket(self, struct):
        """
        Given a dotbracket_string of the same length as the sequence,
        update it with '-' for missing residues and '&' for breakpoints.

        :returns: A string
        """
        return self.parent._missing_into_db(struct)

    def export_missing(self):
        """
        Return the missing residues in a format which can be
        passed to the constructor of a new Sequence object.
        """
        return self.parent._export_missing()

    @property
    def with_modifications(self):
        return _MODIndexer(self)

    def define_length(self, d):
        val = 0
        for i in range(0, len(d), 2):
            if d[i+1]>=d[i]:
                val += sum(1 for _ in self.iter_resids(
                    self.to_resid(d[i]), self.to_resid(d[i + 1])))
                log.debug("Define length of %s with missing incremented to %s", d, val)
        log.debug("Define length of %s with missing is finally %s", d, val)
        return val


class _MODIndexer(_IndexHelper):
    flag = "show_modifications"

    @property
    def with_missing(self):
        return _WMIndexer(self)


class SeqidList(SequenceABC):
    def __init__(self, arg):
        self._list = list(arg)
        self._lookup = {resid: i for i, resid in enumerate(self._list)}
        if len(self._lookup) != len(self._list):
            # duplicate seqids
            c = Counter(self._list)
            for k, amount in c.most_common():
                if amount > 1:
                    log.error("Seq_id %s  occurs %s times!", k, amount)
                else:
                    break
            if c.most_common()[0][1] == 1:
                raise ValueError("len %{} ({}) != ({}) len {}".format(
                                            self._lookup, len(self._lookup),
                                            len(self._list), self._list))
            raise ValueError(
                "Duplicate Seq_id encountered: {}".format(c.most_common()[0][0]))

    def __getitem__(self, i):
        return self._list[i]

    def __len__(self):
        return len(self._list)

    def index(self, elem):
        try:
            return self._lookup[elem]
        except KeyError:
            raise ValueError("{} not in list".format(elem))

    def __ne__(self, other):
        if not isinstance(other, SeqidList):
            return NotImplemented
        return not self==other

    def __eq__(self, other):
        if isinstance(other, list):
            return self._list == other
        elif not isinstance(other, SeqidList):
            return NotImplemented
        return self._list == other._list

    def __hash__(self):
        return hash(self._list)

    def __repr__(self):
        return "SeqidList({})".format(repr(self._list))

class Sequence(object):
    """
    A class holding the RNA sequence.

    It supports 2 indexing conventions (see below), missing residues
    and modified nucleotids.

    The two indexing conventions are:

    *  1-based indexing into the sequence for which structure information is available
       This uses indices of type integer.
       Integer based indexing can not address missing residues!
    *  Indexing using PDB-style numbers (They can start at any positive or negative number and
       may contain insertion codes or gaps.)
       This uses indices of type fgr.RESID
       PDB-Style indices may address missing residues (Reidues for which only sequence
       but no structure information is present).
    """

    def __init__(self, seq, seqids, missing_residues=None, modifications=None):
        """
        :param seq: A string
        :param missing_residues: A list of dictionaries with the following keys:
                "res_name", "chain", "ssseq", "insertion"
        """
        log.debug("Sequence initialized with %s, %s, %s", seq, list(
            map(fgr.resid_to_str, seqids)), missing_residues)
        # Uses 0-based indexing
        i = 0
        self._breaks_after = []
        remaining = seq
        while remaining:
            c, remaining = remaining[0], remaining[1:]
            if c == '&':
                self._breaks_after.append(i - 1)
            else:
                i += 1
        log.debug("Break-points for seq %s are: %s", seq, self._breaks_after)
        self._seq = seq.replace('&', '')
        self._seqids = SeqidList(seqids)
        self._missing_residues = defaultdict(list)
        self._missing_nts = {}
        if missing_residues:
            self._set_missing_residues(missing_residues)  # A dict seq_id:nt
        # resid : modified res_name
        self._modifications = {}
        if modifications:
            self._modifications.update(modifications)

    def _set_missing_residues(self, missing_residues):
        mr, mnts = _sorted_missing_residues(missing_residues)
        log.debug("Setting missing residues to: %s, %s", mr, mnts)
        self._missing_residues = mr
        self._missing_nts = mnts  # A dict seq_id:nt

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
        return [bp + 1 for bp in self._breaks_after]

    def is_valid(self):
        wrong_chars = set(self._seq) - set("AUGCaugc")
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

    def __hash__(self):
        raise AttributeError()

    def __eq__(self, other):
        log.debug("{} =?= {}".format(repr(self), repr(other)))
        log.debug("type of other=%s, type of string literal in this module is %s", type(
            other).__name__, type("").__name__)
        if isinstance(other, type(self)):
            if self._seq != other._seq:
                log.debug("Sequence different!")
                return False
            if self._breaks_after != other._breaks_after:
                log.debug("Breakpoints different!")
                return False
            if self._seqids != other._seqids:
                log.debug("seq_ids different: %s != %s",
                          self._seqids, other._seqids)
                return False
            if self._missing_nts != other._missing_nts:
                if log.isEnabledFor(logging.DEBUG):
                    for key, val in self._missing_nts.items():
                        if key in other._missing_nts:
                            if self._missing_nts[key] != other._missing_nts[key]:
                                log.debug("nt changed for %s from %s to %s", key,
                                          self._missing_nts[key], other._missing_nts[key])
                        else:
                            log.debug("%s missing in other", key)
                    for key in other._missing_nts:
                        if key not in self._missing_nts:
                            log.debug("%s extra in other", key)
                return False
            log.debug("They are equal")
            return True
        return str(self) == other

    def _getitem(self, key, include_missing=False, show_modifications=False):
        log.debug("_getitem called for %s, include_missing=%s, show_modifications=%s",
                  key, include_missing, show_modifications)
        if isinstance(key, int):
            key = to_0_based(key)
            if show_modifications and self._seqids[key] in self._modifications:
                return self._modifications[self._seqids[key]]
            else:
                return self._seq[key]
        elif isinstance(key, fgr.RESID):
            try:
                i = self._seqids.index(key)
            except ValueError:
                if key in self._missing_nts:
                    nt = self._missing_nts[key]
                    if include_missing == False:
                        raise IndexError("No structure available for nucleotide '{}'."
                                         "For look-up including missing residues, use"
                                         "`.with_missing[key]`".format(key))
                    if show_modifications and key in self._modifications:
                        return self._modifications[key]
                    return nt
                error = IndexError(
                    "Nucleotide {} is not part of this RNA".format(key))
                with log_to_exception(log, error):
                    log.error("self._missing_nts = %s", self._missing_nts)
                raise error
            else:
                if show_modifications and key in self._modifications:
                    return self._modifications[key]
                return self._seq[i]
        elif isinstance(key, slice):
            return self._getslice(key, include_missing, show_modifications)
        else:
            raise TypeError("Wrong index type: {}".format(type(key).__name__))

    def _getslice(self, key, include_missing=False, show_modifications=False):
        """
        Do type checks and delegate to self._resid_slice or self._integer_slice
        """
        if key.step is not None and key.step not in [1, -1]:
            raise IndexError("The step in the slice has to be 1 or -1")
        if not isinstance(key.start, (type(None), int, fgr.RESID)):
            raise TypeError(
                "Slice-start must be either an integer, a RESID or None")
        if not isinstance(key.stop, (type(None), int, fgr.RESID)):
            raise TypeError(
                "Slice-stop must be either an integer, a RESID or None")
        if isinstance(key.start, fgr.RESID) or (key.start is None and isinstance(key.stop, fgr.RESID)):
            if not isinstance(key.stop, (fgr.RESID, type(None))):
                raise TypeError(
                    "Mixing integers and RESID for slicing is not allowed")
            return self._resid_slice(key, include_missing, show_modifications)
        return self._integer_slice(to_0_based(key), include_missing, show_modifications)

    def _integer_slice(self, key, include_missing=False, show_modifications=False):
        """
        This provides an implementation for include_missing=False
        and delegates to resid_slice otherwise.

        :param key: Already 0-based, no None in slice
        """
        if ((include_missing == False or not self._missing_nts) and
                (show_modifications == False or not self._modifications)):
            seq = self._seq[key]
            if key.step == -1:
                start = key.stop
                if start is not None:
                    start += 1
            else:
                start = key.start
            return _insert_breakpoints_simple(seq, self._breaks_after, start, key.step == -1)
        if key.start is None:
            startRES = None
        else:
            startRES = self._seqids[key.start]
        if key.stop is None:
            stopRES = None
        else:
            if key.step == -1:
                stopRES = self._seqids[key.stop + 1]
            else:
                if key.stop > 0:
                    stopRES = self._seqids[key.stop - 1]
                else:
                    stopRES = self._seqids[key.stop]
        key2 = slice(startRES, stopRES, key.step)
        log.debug("Converted integer_slice %s to resid-slice %s", key, key2)
        return self._resid_slice(key2, include_missing, show_modifications)

    def _resid_slice_to_int_slice(self, key):
        if key.step == -1:
            default_start = len(self._seqids)
            default_stop = None
        else:
            default_start = 0
            default_stop = len(self._seqids)
        start_i, _ = self._resid_to_index(key.start, default_start, False)
        stop_i, _ = self._resid_to_index(key.stop, default_stop, False)
        if key.step == -1:
            if stop_i is not None:
                stop_i -= 1
        else:
            stop_i += 1
        key = slice(start_i, stop_i, key.step)
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

    def _resid_slice(self, key, include_missing, show_modifications):
        """
        This provides an implementation for inclue_missing=True
        and delegates to integer_slice otherwise.

        :param key: A validated resid slice.
        """
        if include_missing or show_modifications:
            return self._resid_slice_main(key, include_missing, show_modifications)
        else:
            return self._resid_slice_to_int_slice(key)

    def _resid_slice_main(self, key, with_missing, show_modifications):
        if key.step == -1:
            start, stop = key.stop, key.start
        else:
            start, stop = key.start, key.stop
        seqids = list(self.iter_resids(start, stop, with_missing))
        if key.step == -1:
            seqids.reverse()
        seqs = [[]]
        seq = seqs[-1]  # a reference
        for seqid in seqids:
            if seqid == "&":
                if show_modifications:
                    seqs.append([])
                    seq = seqs[-1]
                else:
                    seq.append("&")
            else:
                if show_modifications and seqid in self._modifications:
                    seq.append(self._modifications[seqid])
                else:
                    try:
                        seq.append(self._seq[self._seqids.index(seqid)])
                    except ValueError:
                        seq.append(self._missing_nts[seqid])
        if not show_modifications:
            assert seqs == [seq]
            seqs = "".join(seq)
        return seqs

    def _missing_into_db(self, dotbracket):
        db = ""
        check_i = 0
        for seqid in self.iter_resids(None, None, True):
            if seqid == "&":
                db += "&"
            else:
                try:
                    i = self._seqids.index(seqid)
                    assert i == check_i
                except ValueError:
                    db += "-"
                else:
                    db += dotbracket[i]
                    check_i += 1
        return "".join(db)

    def iter_resids(self, start, stop, include_missing):
        left_res = None
        start_i, start_is_missing = self._resid_to_index(start, 0, True)
        stop_i, stop_is_missing = self._resid_to_index(
            stop, len(self._seqids), True)
        if (start_is_missing or stop_is_missing) and not include_missing:
            raise IndexError(
                "Start or stop are missing residues, but indexing without missing residues was requested.")

        if not start_is_missing and start is not None:
            left_res = start

        # Flag used later if missing_residues == True
        found_start = not start_is_missing
        seq = ""
        if left_res is not None:
            yield left_res
            # seq+=self._seq[self._seqids.index(left_res)]
            #log.debug("Added first residue %s: %s", left_res, seq)
            start_i += 1
        log.debug("Iterating over sequence %s-%s", start_i, stop_i + 1)
        mr_keys = None
        for i in range(start_i, stop_i + 1):
            try:
                right_res = self._seqids[i]
            except IndexError:
                right_res = None
            # Some consistency checks:
            if left_res is not None and right_res is not None and left_res.chain != right_res.chain:
                if i - 1 not in self._breaks_after:
                    raise ValueError("Inconsist break-points discovered."
                                     "Missing breakpoint at position {} "
                                     "between chain {} and chain {}. "
                                     "Break points are {}".format(i - 1,
                                                                  left_res.chain, right_res.chain, self._breaks_after))
            if i - 1 in self._breaks_after:
                if left_res is not None and right_res is not None and left_res.chain == right_res.chain:
                    raise ValueError("Inconsist break-points discovered."
                                     "Breakpoint at position {} "
                                     "in the middle of chain {}".format(i - 1,
                                                                        left_res.chain))
            if include_missing:
                mr, mr_keys = self._missing_residues_between(
                    left_res, right_res)
                old_r = left_res
                for j, r in enumerate(mr_keys):
                    if not found_start and r == start:
                        log.debug("Found start for i=%d, %s", i, r)
                        found_start = True
                        old_r = None
                    log.debug("Now processing missing residue %s", r)
                    if found_start:
                        if old_r is not None and old_r.chain != r.chain:
                            yield "&"
                            log.debug("Breakpoint inserted!")
                        log.debug("Adding missing residue %s", mr[j])
                        yield r
                    log.debug("Comparing %s and stop=%s", r, stop)
                    if r == stop:
                        log.debug("Stop found! Seq=%s", seq)
                        return
                    old_r = r
            # If key.start and key.stop resp. are not modified residues,
            # start_i and stop_i are already accurrate, making
            # additional checks here unneccessary.
            if found_start and i < len(self._seq):
                log.debug("Adding regular residue %s", self._seq[i])
                if mr_keys:  # No missing residues
                    if r.chain != right_res.chain:
                        yield "&"
                        log.debug("Breakpoint inserted!")
                else:
                    if left_res is not None and right_res is not None and left_res.chain != right_res.chain:
                        yield "&"
                        log.debug("Breakpoint inserted!")
                yield right_res
            left_res = right_res
            log.debug("seq now is %s", seq)
        # If start or stop are not part of sequence,
        # make sure they were seen in missing residues
        if not found_start or stop_is_missing:
            raise IndexError("At least one of the residues {} and {} "
                             "is not part of this RNA".format(start, stop))
        return

    def _missing_residues_between(self, from_resid, to_resid):
        log.debug("Searching missing residues between %s  and %s",
                  from_resid, to_resid)
        if from_resid is not None and to_resid is not None:
            if from_resid.chain != to_resid.chain:
                left = self._missing_residues_between(from_resid, None)
                right = self._missing_residues_between(None, to_resid)
                return left[0] + right[0], left[1] + right[1]
        if from_resid is not None:
            from_key = _resid_key(from_resid)
            chain = from_resid.chain
        else:
            from_key = _Smallest()
        if to_resid is not None:
            to_key = _resid_key(to_resid)
            chain = to_resid.chain
        else:
            to_key = _Biggest()
        seq = ""
        res_ids = []
        for res1 in self._missing_residues[chain]:
            if _resid_key(res1) > from_key and _resid_key(res1) < to_key:
                seq += self._missing_nts[res1]
                res_ids.append(res1)
        return seq, res_ids

    def to_resid(self, i):
        return self._seqids[i - 1]

    def to_integer(self, resid):
        return self._seqids.index(resid) + 1

    def iter_modifications(self):
        """
        Iter over tuples seq_id, modified residue code
        """
        for k in sorted(self._modifications):
            yield k, self._modifications[k]

    @property
    def with_missing(self):
        """
        Allows slices with missing residues!
        """
        if self._seqids:
            return _WMIndexer(self)
        else:  # For performance reasons
            return self

    @property
    def with_modifications(self):
        if self._modifications:
            return _MODIndexer(self)
        else:
            return self

    def get_bg_str(self):
        """
        Used during bg-file creation
        """
        out = []
        out.append("seq {}".format(str(self)))
        out.append("seq_ids {}".format(
            " ".join(map(fgr.resid_to_str, self._seqids))))
        for resid, nt in self._missing_nts.items():
            out.append(MissingResidue(resid, nt).to_bg_string())
        for resid, label in self._modifications.items():
            out.append("modification {} {}".format(
                fgr.resid_to_str(resid), label))
        return "\n".join(out) + "\n"

    def define_length(self, d):
        val = 0
        for i in range(0, len(d), 2):
            val += d[i + 1] - d[i] + 1
        log.debug("Define length of %s without missing is %s", d, val)

        return val

    def _export_missing(self):
        """
        Implementation for Sequence().with_missing.export_missing
        """
        out = []
        for resid, nt in self._missing_nts:
            out.append({"RESID": resid, "res_name": nt})
        return out

    def __add__(self, other):
        #print("Sequence add called")
        if isinstance(other, Sequence):
            raise NotImplementedError(
                "This is not yet implemented (but really should be)")
        elif isinstance(other, str):
            return str(self) + other
        else:
            return NotImplemented

    def __radd__(self, other):
        #print("Sequence radd called, other is ", repr(other), type(other).__name__, str.__name__)
        #print(isinstance(other, str))
        if isinstance(other, Sequence):
            raise NotImplementedError(
                "This is not yet implemented (but really should be)")
        elif isinstance(other, str):
            #print("Other is str")
            return other + str(self)
        else:
            return NotImplemented


class SequenceLoader:
    """
    Load the sequence-related part during loading of bg-files
    """

    def __init__(self):
        self.mod = {}
        self.mr = []
        self.seq = None
        self.seq_ids = []

    def consume_fields(self, parts):
        """
        Read one splitted line from the forgi file.

        Returns True, if the line was used/ understood, False otherwise
        """
        mr = MissingResidue.from_bg_fields(parts)
        if mr is not None:
            self.mr.append(mr)
            return True
        elif parts[0] == "seq":
            if self.seq is not None:
                raise ValueError("More than one seq-line encountered.")
            self.seq = parts[1]
            return True
        elif parts[0] == "seq_ids":
            if self.seq_ids:
                raise ValueError("More than one seq-ids line encountered.")
            self.seq_ids = list(map(fgr.resid_from_str, parts[1:]))
            return True
        elif parts[0] == "modification":
            self.mod[fgr.resid_from_str(parts[1])] = " ".join(parts[2:])
            return True
        return False

    @property
    def sequence(self):
        if self.seq is None and not self.seq_ids:
            # The breakpoints are only persisted at the seq and seq_id level.
            raise ValueError("Parsing incomplete: No seq found.")
        elif self.seq is None:
            old_chain = None
            seq = ""
            for resid in self.seq_ids:
                if old_chain is not None and resid.chain != old_chain:
                    seq += "&"
                seq += "N"
                old_chain = resid.chain
            self.seq = seq
        elif not self.seq_ids:
            self.seq_ids = _seq_ids_from_seq_str(self.seq)
        return Sequence(self.seq, self.seq_ids, self.mr, self.mod)


def _seq_ids_from_seq_str(seq):
    """
    Get a list of seq_ids with the same length as the sequence,
    respecting cutpoints.

    We start with seq_id A:1. If a '&' is encountered in the sequence,
    a new chain-letter is used.

    :param seq: A string, optionally containing '&' characters.
    """
    seq_strs = seq.split('&')
    seq_ids = []
    for i, seq_str in enumerate(seq_strs):
        for j, s in enumerate(seq_str):
            seq_ids += [fgr.resid_from_str(
                "{}:{}".format(VALID_CHAINIDS[i], j + 1))]
    return seq_ids


def missing_to_normal(seq):
    """
    Return a Sequence object which has all missing residues moved to normal residues.

    :param seq: A forgi.graph.sequence.Sequence object
    """
    new_seq = seq.with_missing[:]
    new_seqids = list(seq.with_missing.iter_resids(None,None))
    return Sequence(new_seq, new_seqids, None, seq._modifications)
