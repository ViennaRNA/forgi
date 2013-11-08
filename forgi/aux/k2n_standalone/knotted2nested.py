#!/usr/bin/env python
# knotted2nested.py

"""Script to remove pseudoknots from pseudoknotted structures.
Pseudoknot removal is the remaval of one or more base pairs from
a structural model (=list of base pairs), such that all the remaining
base pairs in the structure are nested.

All procedures are described in:
S. Smit, K. Rother, J. Heringa, and R. Knight
Manuscript in preparation.

FORMATS (-f and -F parameter)
The two input formats allowed are bpseq and ct. 
    BPSEQ format (as on the Comparative RNA website):
    4 header lines
    sequence_position residue partner_position 
    NOTE: If base is unpaired, partner is 0

    CT format:
    sequence_position residue seq_pos-1 seq_pos+1 partner_pos seq_pos
    NOTE: If base is unpaired, partner is 0

We support also VIENNA format as output. Standard dot-bracket string. Same
length as the corresponding sequence.

METHODS (-m parameter):
EC -- conflict elimination, most conflicting region gets removed first
EG -- conflict elimination, region with worst gain get removed first
IO -- incremental by order
IL -- incremental by length
IR -- incremental by range
OA -- all optimal solutions depending on optimization function and goal
OSR -- single optimal solution picked at random, also depend on function/goal
OSP -- single optimal solution picked by scoring certain properties, also
    depends on function/goal
NR -- Nussinov algorithm restricted to given list of base pairs

OPTIMIZATION FUNCTIONS (-o parameter)
BPS -- maximize the number of base pairs
HB -- maximize the number of hydrogen bonds in Watson-Crick and Wobble pairs
    (GC pair scores 3, AU and GU pairs score 2)

Author: Sandra Smit (S.Smit@few.vu.nl)

Revision History:
File created on 20 Sep 2007.

"""
from __future__ import division
from sys import exit
from rna2d import Pairs
from optparse import OptionParser
from bpseq import BpseqParser, BpseqParseError
from ct_simple import MinimalCtParser, CtError
from output import ct_output, bpseq_output, vienna_output
from knots import conflict_elimination, ec, eg,\
    inc_order, inc_length, inc_range,\
    opt_all, opt_single_random, opt_single_property,\
    num_bps, hydrogen_bonds, nussinov_restricted

KNOWN_INPUT_FORMATS = {'bpseq': BpseqParser, 'ct': MinimalCtParser}
KNOWN_OUTPUT_FORMATS = {'bpseq': bpseq_output, 'ct': ct_output,
    'vienna': vienna_output}
KNOWN_METHODS = {'EC': ec, 'EG': eg, 
    'IO': inc_order, 'IL': inc_length, 'IR': inc_range,\
     'OA': opt_all, 'OSR': opt_single_random, 'OSP': opt_single_property,\
    'NR': nussinov_restricted,}
KNOWN_OPT_METHODS = {'BPS': ('max', num_bps),\
    'HB': ('max', hydrogen_bonds)}

DEFAULT_INPUT_FORMAT = 'ct'
DEFAULT_OUTPUT_FORMAT = None
DEFAULT_METHOD = 'EG'
DEFAULT_OPT_METHOD = 'BPS'
DEFAULT_VERBOSE = False
DEFAULT_REMOVED = False

class KnottedStructure(object):
    """Lightweight object to hold a pseudoknotted structure
    """

    def __init__(self, Knotted, Seq=None, Header=None):
        """Initialize new KnottedStructure object

        Knotted -- Pairs object representing a (knotted) structure
        Seq -- RNA sequence
        Header -- list of header lines
        """
        self.Knotted = Knotted
        self.Seq = Seq
        self.Header = Header

class NestedStructure(object):
    """Lightweight object to hold a pseudoknot-free structure
    """

    def __init__(self, Nested, Seq=None, Header=None, Removed=None):
        """Initialize new NestedStructure object

        Nested -- Pairs object representing a nested structure
        Seq -- RNA sequence
        Header -- list of header lines
        Removed -- Pairs object containing the removed base pairs
        """
        self.Nested = Nested
        self.Seq = Seq
        self.Header = Header
        self.Removed = Removed


def parse_command_line_parameters():
    """Parse command line arguments"""
    
    # USAGE AND HELP INFO
    usage = '\n'.join(["Usage: python %prog [options] structure_file",\
        "See web or script documentation for an explanation"+\
        " on the formats and methods."])
    input_format_help =\
        "Input format. Options: [%s]. DEFAULT = ct."\
        %(', '.join(KNOWN_INPUT_FORMATS.keys()))
    output_format_help =\
        "Output format. Options: [%s]. DEFAULT = same as input format."\
        %(', '.join(KNOWN_OUTPUT_FORMATS.keys()))
    method_help = "Pseudoknot removal method. Options: [%s]. DEFAULT = EG."\
        %(', '.join(KNOWN_METHODS.keys()))
    opt_method_help="Optimization method. Options: [%s]. DEFAULT = BPS."\
        %(', '.join(KNOWN_OPT_METHODS.keys()))

    result = OptionParser(usage=usage)

    # FLAGS
    result.add_option('-v','--verbose',action='store_true',dest='verbose',\
        help="Run script with verbose output")
    result.add_option('-r','--removed',action='store_true',dest='removed',\
        help="Report removed base pairs")
    # VLAUE PARAMETERS
    result.add_option('-f','--input_format',action='store',type='string',\
        dest='input_format', metavar="FORMAT",\
        help=input_format_help)
    result.add_option('-F','--output_format',action='store',type='string',\
        dest='output_format', metavar="FORMAT",\
        help=output_format_help)
    result.add_option('-m','--method', action='store', type='string',\
        dest='method', metavar='METHOD',\
        help=method_help)
    result.add_option('-o','--opt_method', action='store', type='string',\
        dest='opt_method', help=opt_method_help)

    # DEFAULT VALUES 
    result.set_defaults(verbose=DEFAULT_VERBOSE)
    result.set_defaults(removed=DEFAULT_REMOVED)
    result.set_defaults(input_format=DEFAULT_INPUT_FORMAT)
    result.set_defaults(output_format=DEFAULT_OUTPUT_FORMAT)
    result.set_defaults(method=DEFAULT_METHOD)
    result.set_defaults(opt_method=DEFAULT_OPT_METHOD)
    
    return result.parse_args()
    

def k2n_main(input_filename, input_format=DEFAULT_INPUT_FORMAT,\
    output_format=DEFAULT_OUTPUT_FORMAT, method=DEFAULT_METHOD,\
    opt_method=DEFAULT_OPT_METHOD, verbose=DEFAULT_VERBOSE,\
    removed=DEFAULT_REMOVED):
    """Return string of output containing nested structures (+commens)
    
    Use help function or see script documentation for an explanation
        of the parameters.
    """
    result = [] #list of lines with output

    # by default, output is in same format as input
    if output_format is None:
        output_format = input_format

    # =========================================================================
    # VALIDATE INPUT INFORMATION
    # =========================================================================
    # validate input, give feedback on invalid input
    if input_format not in KNOWN_INPUT_FORMATS:
        message = "'%s' is an invalid format. Valid options are [%s]"\
            %(format, ', '.join(KNOWN_INPUT_FORMATS.keys()))
        raise ValueError(message)
    if output_format not in KNOWN_OUTPUT_FORMATS:
    	message = "'%s' is an invalid output format. Valid options are [%s]"\
            %(format, ', '.join(KNOWN_OUTPUT_FORMATS.keys()))
        raise ValueError(message)
    if method not in KNOWN_METHODS:
        message = "'%s' is an invalid method. Valid options are [%s]"\
            %(method, ', '.join(KNOWN_METHODS.keys()))
        raise ValueError(message)
    if opt_method not in KNOWN_OPT_METHODS:
        message =\
            "'%s' is an invalid optimization function. Valid options are [%s]"\
            %(opt_method, ', '.join(KNOWN_OPT_METHODS.keys()))
        raise ValueError(message)
    
    if verbose:
        result.append('INPUT FILE = %s;'%(input_file))
        result.append('FORMAT = %s; METHOD = %s;'%(input_format, method))

    # =========================================================================
    # COLLECT FILE INFORMATION
    # =========================================================================
    # For bpseq files, put info in a list, because it is always only one record
    input_parser = KNOWN_INPUT_FORMATS[input_format]
    input_data = []

    try:
        if input_format == 'bpseq':
            header, seq, pairs = input_parser(open(input_filename,'U'))
            ks = KnottedStructure(pairs, Seq=seq, Header=header)
            input_data = [ks]
        else:
            for header, seq, pairs in input_parser(open(input_filename,'U')):
                input_data.append(KnottedStructure(Knotted=pairs, Seq=seq,\
                    Header=header))
    except (BpseqParseError, CtError, ValueError), e:
        message_lines = ["Cannot read input file.",
            "Make sure the input file is in the specified format.",
            "Internal error message:","%s."%(str(e))]
        raise ValueError('\n'.join(message_lines))
    # =========================================================================
    # CREATE NESTED STRUCTURE(S)
    # =========================================================================
    if not input_data:
        raise ValueError("No input data found")

    pk_function = KNOWN_METHODS[method]
    
    output_data = []
    for knotted_struct in input_data:
        if method in ['OA', 'OSR','OSP']:
            if opt_method == 'HB':
                opt_goal, opt_wrapper =\
                    KNOWN_OPT_METHODS[opt_method]
                opt_function = opt_wrapper(seq)
            else:
                opt_goal, opt_function =\
                    KNOWN_OPT_METHODS[opt_method]
            if verbose:
                result.append('OPTIMIZATION GOAL = %s; FUNCTION = %s;'\
                    %(opt_method, opt_goal))
        else:
            if verbose:
                result.append('OPTIMIZATION GOAL AND FUNCTION IGNORED')
        
        if method == 'OA': # return value is a list of Pairs objects
            nested_result = pk_function(knotted_struct.Knotted,\
                goal=opt_goal, scoring_function=opt_function,\
                return_removed=removed)
            if removed:
                for nested, removed_pairs in nested_result:
                    ns = NestedStructure(Nested=nested, Seq=knotted_struct.Seq,\
                    Header=knotted_struct.Header, Removed=Pairs(removed_pairs))
                    output_data.append(ns)
                    if verbose:
                        result.append(\
                        'NUMBER OF BASE PAIRS BEFORE: %s; AFTER: %s'\
                        %(len(knotted_struct.Knotted), len(ns.Nested)))
            else:
                for nested in nested_result:
                    ns = NestedStructure(Nested=nested, Seq=knotted_struct.Seq,\
                        Header=knotted_struct.Header)
                    output_data.append(ns)
                    if verbose:
                        result.append(\
                        'NUMBER OF BASE PAIRS BEFORE: %s; AFTER: %s'\
                        %(len(knotted_struct.Knotted), len(ns.Nested)))
        elif method in ['OSR','OSP']:
            if removed:
                nested_pairs, removed_pairs = pk_function(pairs, goal=opt_goal,\
                    scoring_function=opt_function, return_removed=removed)
                ns = NestedStructure(Nested=nested_pairs,\
                    Seq=knotted_struct.Seq,\
                    Header=knotted_struct.Header, Removed=Pairs(removed_pairs))
            else:
                nested_pairs = pk_function(pairs, goal=opt_goal,\
                    scoring_function=opt_function)
                ns = NestedStructure(Nested=nested_pairs,\
                    Seq=knotted_struct.Seq,\
                    Header=knotted_struct.Header)
            if verbose:
                result.append('NUMBER OF BASE PAIRS BEFORE: %s; AFTER: %s'\
                    %(len(pairs), len(nested_pairs)))
            output_data.append(ns)
        else: # return value is a single Pairs object
            if removed:
                nested_pairs, removed_pairs = pk_function(pairs,\
                    return_removed=removed)
                ns = NestedStructure(Nested=nested_pairs,\
                    Seq=knotted_struct.Seq,\
                    Header=knotted_struct.Header, Removed=Pairs(removed_pairs))
            else:
                nested_pairs = pk_function(pairs)
                ns = NestedStructure(Nested=nested_pairs,\
                    Seq=knotted_struct.Seq,\
                    Header=knotted_struct.Header)
            if verbose:
                result.append('NUMBER OF BASE PAIRS BEFORE: %s; AFTER: %s'\
                    %(len(pairs), len(nested_pairs)))
            output_data.append(ns)

    # =========================================================================
    # PRODUCE OUTPUT
    # =========================================================================
    
    output_function = KNOWN_OUTPUT_FORMATS[output_format]
    
    # return string with output
    for ns in output_data:
        result.append(output_function(seq=ns.Seq, pairs=ns.Nested,\
            header_lines=ns.Header, removed=ns.Removed))
    return '\n'.join(result)

if __name__ == "__main__":
     
    # =========================================================================
    # PARSE COMMAND LINE PARAMETERS AND SEND INFO TO MAIN FUNCTION
    # =========================================================================
    
    opts,arg = parse_command_line_parameters()

    try:
        input_file = arg[0]
    except IndexError:
        raise ValueError("No input file specified")
            
    result = k2n_main(input_file, input_format=opts.input_format,\
    output_format=opts.output_format, method=opts.method,\
    opt_method=opts.opt_method, verbose=opts.verbose, removed=opts.removed)
    print result
