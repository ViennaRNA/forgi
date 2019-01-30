import copy
import re
import forgi.graph.bulge_graph as fgb
import forgi.graph.sequence as fgs
import forgi.graph.residue as fgr
import forgi.threedee.model.coarse_grain as ftmc
import logging
import traceback
import numpy as np

log = logging.getLogger(__name__)


class NumberedDotbracket():
    def __init__(self, dotbracket_str="", residue_numbers=None, helix_ends=None):
        """
        :param dotbracket_str: A string
        :param residue numbers: A list of integers or forgi RESIDS.
                                Same length as dotbracket_str
        :param helix_ends: If a basepair represents a helix or
                           stem with missing residues,
                           The helix end give the last/ first missing(condensed)
                           residue number which is represented by this basepair.
                           I.e. the basepair of the helix that is at the 3' end
                           of the forward strand.
                           It is a list with the same length as residue_numers.

                           *   Example 1. Every nucleotide only represents itself::
                                    dotbracket_str   (((..)))
                                    residue numbers  12345678
                                    helix_ends       12345678
                           *  Example 2. The same helix as Ex. 1, now condensed.
                                    The helix 1-3,6-8 is represented by the basepair 1-8::
                                    dotbracket_str   (..)
                                    residue numbers  1458
                                    helix_ends       3456
                                    Here the nucleotide 1 represents one strand
                                    of the helix, going from 1 to 3.
                                    The nucleotide 8 represents the other strand,
                                    going from 6 to 8.

                          The value for helix ends for unpaired nucleotides is
                          undefined and should ignored.
        """
        dotbracket_str = dotbracket_str.replace("&","")
        if residue_numbers is None:
            residue_numbers=[]
        assert len(dotbracket_str)==len(residue_numbers)
        self.dotbracket_str = dotbracket_str
        self.residue_numbers = residue_numbers
        if helix_ends == None:
            self.helix_ends = copy.copy(self.residue_numbers)
        else:
            self.helix_ends = helix_ends
        if False: #self.residue_numbers:
            tb = traceback.format_list(traceback.extract_stack()[:-2])
            for line in tb:
                log.debug(line.strip())
            log.debug("NumberedDotbracket initialized. "
                  "residue_numbers[0]=%s", residue_numbers[0])

    def without_unpaired(self):
        """
        Return NumberedDotbracket with dots removed

        :param dotbracket: A dotbracket string
        :param residue_numbers: A list of residue numbers,
                                with same length as dotbracket
        """
        log.debug("Removing unpaired from %s",
                  list(map(fgr.resid_to_str, self.residue_numbers)))
        new_dotbracket = ""
        new_residue_numbers = []
        new_helix_ends = []

        for i, ch in enumerate(self.dotbracket_str):
            if ch!='.':
                new_dotbracket+=ch
                new_residue_numbers.append(self.residue_numbers[i])
                new_helix_ends.append(self.helix_ends[i])
        log.debug("Residue numbers now %s",
                  list(map(fgr.resid_to_str, new_residue_numbers)))
        return NumberedDotbracket(new_dotbracket, new_residue_numbers, new_helix_ends)


    def _to_forgi(self):
        try: # Expects residue_numers to be forgi.graph.residue.RESID instances
            breakpoints = ftmc.breakpoints_from_seqids(self.residue_numbers)
        except AttributeError: #residue numbers are ints or something else
            breakpoints=[]
        dotbracketstr = fgs._insert_breakpoints_simple(self.dotbracket_str, breakpoints, 1)
        log.info("Converting %s to forgi", dotbracketstr)
        bg = fgb.BulgeGraph.from_dotbracket(dotbracketstr)
        log.debug("Forgi-db is %s ", bg.to_dotbracket_string())

        return bg

    def without_short_helices(self, min_length):
        """
        Check the stems if the correspond to the minimum length given by the Users

        :param min_length: minimum stem length
        """
        bg = self._to_forgi()
        to_unpaired = []
        for stem in bg.stem_iterator():
            if bg.stem_length(stem)<min_length:
                to_unpaired.extend(bg.define_residue_num_iterator(stem))
        numpy_db = np.array(list(self.dotbracket_str))
        log.info("%s %s", "".join(numpy_db), type(numpy_db))
        to_unpaired = np.array(to_unpaired, dtype = int)
        to_unpaired-=1 #0-based indexing, numpy supports element-wise substraction
        numpy_db[to_unpaired]="."
        log.info("After to_unpaired %s ", "".join(numpy_db))
        return NumberedDotbracket("".join(numpy_db), self.residue_numbers, self.helix_ends)

    def condensed(self, remove_helices_shorter_than=3):
        """
        Return input with consecutive basepairs remove_unpaired

        :param dotbracket: A dotbracket string
        :param residue_numbers: A list of residue numbers,
                                with same length as dotbracket
        """
        log.debug("Condensing: Seq-ids are %s of length %s",
                  list(map(fgr.resid_to_str, self.residue_numbers)),
                  len(self.residue_numbers))
        bg = self._to_forgi()
        log.debug("Breakpoints are (0-based) %s", bg.seq._breaks_after)
        log.debug("Before setting: ResNum: %s will be set to %s",
                  bg.seq._seqids, self.residue_numbers)
        bg.seq._seqids = self.residue_numbers
        log.debug("After setting, before condensing: ResNum: %s", bg.seq._seqids[0])
        bg2 = bg.transformed.condensed()
        log.debug("Successfully transformed")
        log.debug("DB now %s", bg2.to_dotbracket_string())
        log.debug("Resids now %s", bg2.seq._seqids)

        helix_ends = []
        for seqid in  bg2.seq._seqids:
            elem_original = bg.get_elem(seqid)
            if seqid == bg.seq.to_resid(bg.defines[elem_original][0]):
                helix_end = bg.seq.to_resid(bg.defines[elem_original][1])
            elif seqid == bg.seq.to_resid(bg.defines[elem_original][3]):
                helix_end = bg.seq.to_resid(bg.defines[elem_original][2])
            else:
                assert False
            helix_ends.append(helix_end)

        log.debug("Helix ends: %s", helix_ends)
        return NumberedDotbracket(bg2.to_dotbracket_string(), bg2.seq._seqids, helix_ends)

    def __str__(self):
        return str(self.dotbracket_str)

    def __repr__(self):
        return "NumberedDotbracket({}, {})".format(repr(self.dotbracket_str),
                                                   repr(self.residue_numbers))

    def __contains__(self, value):
        return value in self.dotbracket_str

    def __iter__(self):
        for i, bracket in enumerate(self.dotbracket_str):
            yield NumberedDotbracket(bracket, [self.residue_numbers[i]], helix_ends=[self.helix_ends[i]])

    def __add__(self, other):
        db = self.dotbracket_str +  other.dotbracket_str
        rn = self.residue_numbers + other.residue_numbers
        he = self.helix_ends + other.helix_ends
        log.debug("Add to PK: %r, %r", self, other)
        return NumberedDotbracket(db, rn, helix_ends=he)

    def __eq__(self, other):
        if isinstance(other, NumberedDotbracket):
            return (self.dotbracket_str==other.dotbracket_str and
                    self.residue_numbers == other.residue_numbers)
        else:
            # If other is a string, compare it to self.dotbracket_str
            return self.dotbracket_str == other

    def __getitem__(self, key):
        if isinstance(key, int):
            return NumberedDotbracket(self.dotbracket_str[key], [self.residue_numbers[key]], helix_ends=[self.helix_ends[key]])
        elif isinstance(key, slice):
            return NumberedDotbracket(self.dotbracket_str[key], self.residue_numbers[key], helix_ends=self.helix_ends[key])

    def __hash__(self):
        return hash(self.dotbracket_str)

    def __len__(self):
        return len(self.dotbracket_str)

    def count(self, pattern):
        return (self.dotbracket_str.count(pattern))

    def without_substr(self, pattern):
        """
        Return number of occurrences of pattern and
        a NumberedDotbracket all occurrences removed.

        Occurrences of the pattern are removed recursicvely. E.g.:
        The pattern "ab" will be removed from "aabb" two times.
        E.g.: two occurrences of "([)]" will be removed from "(([)][)]"
        """
        db =self.dotbracket_str
        res_nums = copy.copy(self.residue_numbers)
        h_ends = copy.copy(self.helix_ends)
        removed = []
        while True:
            try:
                i = db.index(pattern)
            except ValueError:
                return removed, NumberedDotbracket(db, res_nums)
            removed.append(NumberedDotbracket(db[i:i+len(pattern)],
                                              res_nums[i:i+len(pattern)],
                                              h_ends[i:i+len(pattern)]))
            db = db[:i]+db[i+len(pattern):]
            res_nums = res_nums[:i]+res_nums[i+len(pattern):]
            h_ends = h_ends[:i]+h_ends[i+len(pattern):]
