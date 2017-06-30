#!/usr/bin/python

from __future__ import print_function
import sys, math, os
import forgi.graph.bulge_graph as cgb
import forgi.utilities.debug as cud
import forgi.utilities.stuff as cus

from optparse import OptionParser

def print_rosetta_constraints(bg):
    '''
    Convert the stems into base pair constraints compatible
    with the rna_denovo program.

    @param bg: A BulgeGraph data structure
    '''
    for s in bg.stem_iterator():
        for i in range(bg.stem_length(s)):
            print("STEM PAIR %d %d" % (bg.defines[s][0] + i, bg.defines[s][3] - i))

def main():
    usage = """
        Usage: ./graph_to_rosetta_adjacent_stems.py

        Take each set of stems that are separated by something and create
        a set of rosetta constraings, such that they can be folded keeping 
        the intermediate structure intact.
        """
    parser = OptionParser(usage = usage)
    parser.add_option("-o", "--output-dir", dest="output_dir", default=".", help="Place the resulting sequence and constraint files in this directory", type='str')
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)
    if args[0] == '-':
        f = sys.stdin
    else:
        f = open(args[0])

    bg_str = "".join(f.readlines())

    bg = cgb.BulgeGraph()
    bg.from_bg_string(bg_str)

    for (s1, e1, s2) in bg.adjacent_stem_pairs_iterator():
        seq_str = ""
        param_str = ""

        # group the defines into pairs of intervals
        pairs_s1 = zip(*[iter(bg.defines[s1])]*2)
        pairs_s2 = zip(*[iter(bg.defines[s2])]*2)
        pairs_e1 = zip(*[iter(bg.defines[e1])]*2)
        
        # find the longest contiguous stretches of nucleotides
        intervals = pairs_s1 + pairs_e1 + pairs_s2
        merged_intervals = cus.merge_intervals(intervals)

        #create a lookup table relating the nucleotide numbers in the
        #bulge graph to the nucleotide numbers in the newly formed
        #sequence which will be fed to rosetta
        res_num_map = dict()
        res_num = 1
        prev_cutpoint_open = ""
        for i in merged_intervals:
            param_str += prev_cutpoint_open
            seq_str += bg.seq[i[0]-1:i[1]]
            for j in range(i[0], i[1] + 1):
                res_num_map[j] = res_num
                res_num += 1
            prev_cutpoint_open = "CUTPOINT_OPEN %d\n" % (res_num-1)

        # if the separating element is an interior loop, then we should add
        # a cutpoint that will be closed right after it's first nucleotide,
        # which is actually after the last nucleotide of the stem
        if e1[0] == 'i':
            param_str += "CUTPOINT_CLOSED %d\n" % (res_num_map[max(bg.defines[e1])])

        for (bpl, bpr) in bg.stem_bp_iterator(s1):
            param_str += "STEM PAIR %d %d W W A\n" % (res_num_map[bpl], res_num_map[bpr])
        for (bpl, bpr) in bg.stem_bp_iterator(s2):
            param_str += "STEM PAIR %d %d W W A\n" % (res_num_map[bpl], res_num_map[bpr])

        fn_prefix = os.path.join(options.output_dir, "%s_%s_%s/" % (s1, e1, s2))

        # Create the directories that will store the information from the rosetta run
        try:
            os.makedirs(fn_prefix + "/rosetta_inputs")
        except:
            pass

        try:
            os.makedirs(fn_prefix + "/rosetta_outputs")
        except:
            pass

        try:
            os.makedirs(fn_prefix + "/cluster_output")
        except:
            pass

        try:
            os.makedirs(fn_prefix + "/cluster_error")
        except:
            pass

        try:
            os.makedirs(fn_prefix + "/job_info")
        except:
            pass

        # Output the sequence and constraints in the rosetta_inputs directory
        # The other directories will be populated as rosetta runs
        with open(fn_prefix + "rosetta_inputs/add_seq.fasta", 'w') as f:
            f.write(">" + bg.name + "\n")
            f.write(seq_str.lower() + "\n")
        with open(fn_prefix + "rosetta_inputs/add_constraints.cnt", 'w') as f:
            f.write(param_str)

if __name__=="__main__":
    main()
