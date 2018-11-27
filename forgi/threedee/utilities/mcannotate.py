#!/usr/bin/python
from __future__ import print_function, unicode_literals
from builtins import str
from builtins import range
import sys
import copy
import re
import logging
from collections import defaultdict

from logging_exceptions import log_to_exception, log_exception

import forgi.utilities.debug as fud
import forgi.graph.residue as fgr

log = logging.getLogger(__name__)


def parse_resid(mcann_resid):
    '''
    Read in an MC-Annotate formatted resid and return a tuple
    that can be used as an index into a Bio.PDB.Chain.
    '''
    parts = mcann_resid.split('.')
    if len(parts) == 1:
        return (' ', int(parts[0]), ' ')
    else:
        return (' ', int(parts[0]), parts[1])


def format_resid(pdb_resid):
    '''
    Convert a PDB.Chain.Residue id to an MC-Annotate formatted
    residue identifier.
    '''

    if pdb_resid[2] == ' ':
        return str(pdb_resid[1])
    else:
        return str(pdb_resid[1]) + "." + pdb_resid[2]


def parse_chain_base(chain_base):
    """
    Parse the string identifying a chain and a base in an MC-Annotate generated
    annotation file.

    As an example, the string 'A33' means residue 33 on chain 'A'.

    @param chain_base: The MC-Annotate identification string (i.e. 'A33')
    @return: A pair containing the chain id and residue number (i.e. ('A', 33))
    """
    log.debug("Chain_base is '{}'".format(chain_base))
    if (ord(chain_base[0]) >= ord('A') and ord(chain_base[0]) <= ord('Z')) or (ord(chain_base[0]) >= ord('a') and ord(chain_base[0]) <= ord('z')):
        # normal string
        chain = chain_base[0]
        base = chain_base[1:]
    else:
        # quoted string (i.e. ''33'')
        if chain_base[0] == '\'':
            end_quote_idx = chain_base.find('\'', 1)
            chain = chain_base[1:end_quote_idx]
            base = chain_base[end_quote_idx + 1:]
        else:
            # no chain identifier
            chain = ''
            base = chain_base
    log.debug("base is %s", base)
    return (chain, base)


def parse_base_pair_id(base_pair_id):
    """
    Separate the two chain/base identifiers present in the interaction section of
    an MC-Annotate output file.

    @param base_pair_id: The identifier string for the interacting nucleotides (i.e. 'A33-B45')
    @return: 4-tuple containing of the form (chain1, res1, chain2, res2) i.e. ('A', 33, 'B', '45')
    """
    # A number in single quotes or a letter, followed by a (potentially negative) number and
    # potentiallly by an insertion code.
    residue_pattern = r"(?:'\d'|[A-Za-z])-?\d+(?:\.[A-Za-z])?"

    parts = re.findall(residue_pattern, base_pair_id)
    if len(parts) != 2:
        e = ValueError(
            "Invalid interaction in the MC-Annotate file: %s" % base_pair_id)
        with log_to_exception(log, e):
            log.error("Regex matched the following parts: %s", parts)
        raise e
    if "-".join(parts) != base_pair_id:
        raise ValueError(
            "Invalid interaction in the MC-Annotate file: %s" % base_pair_id)

    log.debug("Parts are '{}'".format(parts))
    (from_chain, from_base) = parse_chain_base(parts[0].strip())
    (to_chain, to_base) = parse_chain_base(parts[1].strip())

    return (from_chain, from_base, to_chain, to_base)


def get_interacting_base_pairs(line):
    """
    Return the identification part of an interaction line of an MC-Annotate output file.

    @param line: The entire interaction line
    @return: The first part, containing the identifiers of the interacting entities (i.e. 'A33-B45')
    """
    line_parts = line.split(' ')
    return parse_base_pair_id(line_parts[0])


def iterate_over_residue_list(mcannotate_lines):
    """
    Generator function allowing the iteration over all of the lines in
    the 'Residue conformations' section of an MC-Annotate output file.

    @param mcannotate_lines: All of the lines of an MC-Annotate output file
    @return: Yield only lines in the 'Residue conformations' section
    """
    residue_conf_line = False
    for line in mcannotate_lines:
        if line.find('Residue conformations') == 0:
            residue_conf_line = True
            continue
        if line.find('Adjacent stacking') == 0:
            residue_conf_line = False
            continue
        if residue_conf_line:
            yield line


def iterate_over_interactions(mcannotate_lines):
    """
    Generator function for the iteration over lines in the 'Base-pairs' section of an
    MC-Annotate output file.

    @param mcannotate_lines: All of the lines of an MC-Annotate output file
    @return: Yield only lines in the 'Base-pairs' section
    """
    base_pair_line = False
    for line in mcannotate_lines:
        if line.find("Base-pairs ---") == 0:
            base_pair_line = True
            continue
        if line.find("Residue conformations") == 0:
            base_pair_line = False
            continue
        if base_pair_line:
            log.debug("A Basepair line is '{}'".format(line))
            try:
                (from_chain, from_base, to_chain,
                 to_base) = get_interacting_base_pairs(line)
            except ValueError as ve:
                log_exception(ve, logging.WARNING, with_stacktrace=False)
                continue

            yield line.strip()


def get_dotplot(lines):
    """docstring for get_dotplot"""
    residues = []
    residue_types = []
    bps = defaultdict(lambda: -1)
    bpseq_str = ""

    for line in iterate_over_residue_list(lines):
        parts = line.split(' ')
        residues.append(parse_chain_base(parts[0]))  # A tuple chain, id
        residue_types += [parts[2]]

    paired = set()
    for line in iterate_over_interactions(lines):
        parts = line.split(' ')
        #bond_type = parts[3]
        # if bond_type.find('Ww/Ww') >= 0 or bond_type.find('Ww/Ws') >= 0 or bond_type.find('Ws/Ww') >= 0:
        if ((line.find('Ww/Ww') >= 0 and (line.find('A-U') >= 0 or
                                          line.find('U-A') >= 0 or
                                          line.find('C-G') >= 0 or
                                          line.find('G-C') >= 0)) or
            (line.find('Ws/Ww') >= 0 and line.find('U-G') >= 0) or
                (line.find('Ww/Ws') >= 0 and line.find('G-U') >= 0)):
            # if bond_type.find('Ww/Ww') >= 0:
            # print line
            chain1, base1, chain2, base2 = parse_base_pair_id(parts[0])
            res1 = (chain1, base1)
            res2 = (chain2, base2)
            if res1 in paired or res2 in paired:
                if log.isEnabledFor(logging.WARNING):
                    if res1 in bps:
                        existing = "{} - {}".format(res1, residues[bps[res1]])
                    else:
                        existing = "{} - {}".format(res2, residues[bps[res2]])
                    log.warning(
                        "Base-triple encountered: Ignoring basepair %s - %s, because basepair %s exists", res1, res2, existing)
                continue

            paired.add(res1)
            paired.add(res2)
            try:
                bps[res1] = residues.index(res2)
                bps[res2] = residues.index(res1)
            except ValueError as e:
                with log_to_exception(log, e):
                    log.error(
                        "bps = %s, residues = %s, res1 = %s, res2 = %s", bps, residues, res1, res2)
                raise

    for i in range(len(residue_types)):
        bpseq_str += "%d %s %s\n" % (i + 1,
                                     residue_types[i], bps[residues[i]] + 1)
    seq_ids = _seqids_from_residue_map(residues)
    return bpseq_str, seq_ids


def _seqids_from_residue_map(residue_map):
    """
    Create the list of seq_ids from the list of MC-Annotate identifiers in the
    residue map.
    """
    seq_ids = []
    for i, r in enumerate(residue_map):
        log.debug("Residue %r", r)
        (from_chain, from_base) = r
        seq_ids += [fgr.RESID(from_chain, parse_resid(from_base))]
    log.debug("Seq_ids are %s", seq_ids)
    return seq_ids
