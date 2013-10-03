#!/usr/bin/env python
# output.py

"""
Author: Sandra Smit (Sandra.Smit@colorado.edu)

Revision History:
File created on 25 Sep 2007.

"""
from __future__ import division
from rna2d import Pairs

def ct_output(seq, pairs, header_lines=None, removed=None):
    """Return sequence and structure information in ct format

    seq -- string of sequence characters, no validation
    pairs -- Pairs object or list of tuples
    header_lines -- list of header lines at the top of the file

    Return value is a formatted string
    """
    pairs = Pairs(pairs)
    result = []
    
    if not header_lines: # None or empty
        header = '%s ENERGY = 0.0'%(len(seq))
    else: # take first line. If ct leave as is, otherwise add energy
        head = header_lines[0]
        if 'ENERGY' in head or 'Energy' in head or 'dG' in head:
            header = head.strip()
        else:
            header = '%s ENERGY = 0.0 %s'%(len(seq), head.strip())
    
    result.append(header)
    partners = pairs.toPartners(len(seq))
    for idx, (seq_symbol, partner_idx) in enumerate(zip(seq, partners)):
        if partner_idx is None:
            partner_idx = -1
        result.append('\t'.join(map(str,[idx+1, seq_symbol, idx,\
            idx+2, partner_idx+1, idx+1])))

    if removed is not None:
        result.append('# REMOVED PAIRS FOR: %s'%(header))
        partners = removed.toPartners(len(seq))
        for idx, (seq_symbol, partner_idx) in enumerate(zip(seq, partners)):
            if partner_idx is None:
                partner_idx = -1
            result.append('\t'.join(map(str,[idx+1, seq_symbol, idx,\
                idx+2, partner_idx+1, idx+1])))

    return '\n'.join(result)

def bpseq_output(seq, pairs, header_lines=None, removed=None):
    """Return sequence and structure information in bpseq format

    seq -- string of sequence characters, no validation
    pairs -- Pairs object or list of tuples
    header_lines -- list of header lines at the top of the file

    Return value is a formatted string
    """
    result = []
    if header_lines:
        result.extend(header_lines)
    partners = pairs.toPartners(len(seq))
    for idx, (seq_symbol, partner_idx) in enumerate(zip(seq, partners)):
        if partner_idx is None:
            partner_idx = -1
        result.append(' '.join(map(str,[idx+1, seq_symbol, partner_idx+1])))
    
    if removed is not None:
        result.append('# REMOVED BASE PAIRS')
        partners = removed.toPartners(len(seq))
        for idx, (seq_symbol, partner_idx) in enumerate(zip(seq, partners)):
            if partner_idx is None:
                partner_idx = -1
            result.append(' '.join(map(str,[idx+1, seq_symbol, partner_idx+1])))

    return '\n'.join(result)

def vienna_output(seq, pairs, header_lines=None, removed=None):
    """Return sequence and structure information in vienna format

    seq -- string of sequence characters, no validation
    pairs -- Pairs object or list of tuples
    header_lines -- list of header lines at the top of the file

    Return value is a formatted string:
    > header information
    sequence
    vienna structure
    """
    result = []
    if not header_lines:
        head = '> Nested structure'
    else: # take first
        head = '> ' + header_lines[0].strip()
    result.append(head)
    result.append(seq)
    result.append(pairs.toVienna(len(seq)))
    if removed is not None:
        result.append('>REMOVED BASE PAIRS FOR: %s'%(head))
        try:
            result.append(removed.toVienna(len(seq)))
        except:
            result.append(str(removed))
    return '\n'.join(result)

if __name__ == "__main__":
    pass
