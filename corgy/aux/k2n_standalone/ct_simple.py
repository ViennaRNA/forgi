#!/usr/bin/env python
# ct.py

"""
Author: Sandra Smit (Sandra.Smit@colorado.edu)

Revision History:
File created on 21 Sept 2007.

"""
from __future__ import division
from rna2d import Pairs

class CtError(Exception):
    pass

def is_ct_line(line):
    """Return True if line is header line
    
    line -- str, single line

    Header line is recognized by one of the following words:
    dG, ENERGY, or Energy
    """
    if 'dG' in line:
        return True
    if 'ENERGY' in line or 'Energy' in line:
        return True
    return False

def ct_record_finder(lines):
    """Yield successive ct records"""
    curr = []
    for line in lines:
        if not line.strip():
            continue
        if is_ct_line(line):
            if curr:
                yield curr
            curr = [line.strip()]
        else:
            curr.append(line.strip())
    if curr:
        yield curr
        curr = []

def ct_head_body(lines):
    """Separate first line (header) from the rest (seq/struct info)

    lines -- list of lines or anything that behaves like it, e.g.
        filestream.
    """
    head = lines[0]
    body = lines[1:]
    return head, body

def ct_parse_header(line):
    """Return simple string of header"""
    return line.strip()

def ct_parse_content(lines):
    """Return tuple of (raw_sequence, Pairs object)
   
    lines -- list of lines or anything behaving like it

    This function parses seq/struct lines of files in .ct format.
    Return value is tuple of sequence and Pairs object
        raw_sequence is simple string of characters
        Pairs object will be zero-based. 

    This method should be universal for .ct files from different sources.
    In general they only differ in their header.
    
    Assumes 1-based, and returns 0-based structure object, corresponding
        to the sequence.
    Not yet: (Parses standard 6-column format.)
    """
    # check assumption (1-based numbering)
    try:
        first_char = int(lines[0].strip().split()[0])
        if first_char != 1:
            raise CtError("First content line starts with %s, should be 1"\
                %(first_char))
    except ValueError:
        raise CtError("First content line doesn't start with an int: %s"\
            %(lines[0].strip()))

    seq = []
    struct = Pairs()
    last_line = None
    for line in lines:
        parts = line.strip().split()
        if len(parts) != 6:
            raise CtError("Number of parts in line is %s, expected 6"\
                %(len(parts)))
        seq.append(parts[1])
        up, down = int(parts[0]), int(parts[4])

        # check file consistency
        if last_line is not None and up != last_line +1:
            raise CtError("Numbers out of order: %s after %s"%(up, last_line))
        if int(parts[2]) != up-1 or int(parts[5]) != up:
            raise CtError("Inconsistency in line: %s"%(line.strip()))
        if int(parts[3]) != up+1 and int(parts[3]) != 0:
            raise CtError("Inconsistency in line: %s"%(line.strip()))
        # if everything is ok, append
        if down == 0:
            struct.append((up, None))
        if down != 0:
            struct.append((up, down))
        last_line = up
    if struct.hasConflicts():
        raise CtError("Stucture has conflicts")
    renumbered = Pairs()
    for x,y in struct.directed():
        renumbered.append((x-1, y-1))
    #struct = adjust_base(struct.directed(), -1)
    seq = ''.join(seq)

    return seq, renumbered

def MinimalCtParser(lines, strict=True, header_parser=ct_parse_header,\
    content_parser=ct_parse_content):
    """Return pure data as is"""
    for rec in ct_record_finder(lines):
        try:
            head, body = ct_head_body(rec)
            header_info = header_parser(head)
            seq, pairs = content_parser(body)
            yield [header_info], seq, pairs
        except CtError, e:
            if strict:
                raise CtError(str(e))
            else:
                continue

def CtParser(lines, strict=True, header_parser=ct_parse_header,\
    content_parser=ct_parse_content):
    """Return RNA sequence object etc.
    """
    for header_info, seq, pairs in MinimalCtParser(lines, strict=strict,\
        header_parser=ct_parse_header, content_parser=ct_parse_content):
        try:
            sequence, pairs = supported[program]['INTEGRATION']\
                (header_info, content_info)
            yield sequence, pairs
        except CtError, e:
            if strict:
                raise CtError(str(e))
            else:
                continue

if __name__ == "__main__":
    from sys import argv
    #for lines in CTRecordFinder(open(argv[1],'U')):
    #    print lines

    for res in MinimalCtParser(open(argv[1],'U')):
        print res
