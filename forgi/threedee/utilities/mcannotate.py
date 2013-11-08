#!/usr/bin/python

import sys
import copy

def parse_chain_base(chain_base):
    """
    Parse the string identifying a chain and a base in an MC-Annotate generated
    annotation file.

    As an example, the string 'A33' means residue 33 on chain 'A'.

    @param chain_base: The MC-Annotate identification string (i.e. 'A33')
    @return: A pair containing the chain id and residue number (i.e. ('A', 33))
    """
    if (ord(chain_base[0]) >= ord('A') and ord(chain_base[0]) <= ord('Z')) or (ord(chain_base[0]) >= ord('a') and ord(chain_base[0]) <= ord('z')):
        # normal string
        chain = chain_base[0]
        base = int(chain_base[1:])
    else:
        # quoted string (i.e. ''33'')
        if chain_base[0] == '\'':
            end_quote_idx = chain_base.find('\'', 1)
            chain = chain_base[1:end_quote_idx]
            base = int(chain_base[end_quote_idx+1:])
        else:
            # no chain identifier
            chain = ''
            base = int(chain_base)

    return (chain, base)

def parse_base_pair_id(base_pair_id):
    """
    Separate the two chain/base identifiers present in the interaction section of
    an MC-Annotate output file.

    @param base_pair_id: The identifier string for the interacting nucleotides (i.e. 'A33-B45')
    @return: 4-tuple containing of the form (chain1, res1, chain2, res2) i.e. ('A', 33, 'B', '45')
    """
    parts = base_pair_id.split('-')
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
            try:
                (from_chain, from_base, to_chain, to_base) =  get_interacting_base_pairs(line)
            except ValueError as ve:
                print >> sys.stderr, "ValueError:", ve
                print >> sys.stderr, "line:", line
                continue

            yield line.strip()

# From:
# http://code.activestate.com/recipes/389639/
class DefaultDict(dict):
    """Dictionary with a default value for unknown keys."""
    def __init__(self, default):
        self.default = default

    def __getitem__(self, key):
        if key in self: 
            return self.get(key)
        else:
            ## Need copy in case self.default is something like []
            return self.setdefault(key, copy.deepcopy(self.default))

    def __copy__(self):
        copy = DefaultDict(self.default)
        copy.update(self)
        return copy

def get_dotplot(lines):
    """docstring for get_dotplot"""
    residues = []
    residue_types = []
    bps = DefaultDict(-1)
    output_str = ""

    for line in iterate_over_residue_list(lines):
        parts = line.split(' ')
        residues += [parts[0]]
        residue_types += [parts[2]]

    paired = set()
    for line in iterate_over_interactions(lines):
        parts = line.split(' ')
        bond_type = parts[3]
        if bond_type.find('Ww/Ww') >= 0 or bond_type.find('Ww/Ws') >= 0 or bond_type.find('Ws/Ww') >= 0:
        #if bond_type.find('Ww/Ww') >= 0:
            parts1 = parts[0].split('-')
            #print line

            if parts1[0] in paired or parts1[1] in paired:
                print >>sys.stderr, "paired:", parts1[0], parts1[1]
                continue

            paired.add(parts1[0])
            paired.add(parts1[1])

            bps[parts1[0]] = residues.index(parts1[1])
            bps[parts1[1]] = residues.index(parts1[0])

    for i in range(len(residue_types)):
        output_str += "%d %s %s\n" % ( i+1, residue_types[i], bps[residues[i]]+1)

    return output_str
