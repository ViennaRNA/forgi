#!/usr/bin/env python
from __future__ import print_function

import logging
import os.path
import sys
import numpy as np

import forgi.graph.bulge_graph as fgb
import forgi.utilities.commandline_utils as fuc

log = logging.getLogger(__name__)


def to_bg_or_cg_string(cg):
    if not hasattr(cg, "coords") or not cg.coords.is_filled:
        return cg.to_bg_string()
    else:
        return cg.to_cg_string()
def cg_or_bg_extention(cg):
    if not hasattr(cg, "coords") or not cg.coords.is_filled:
        return ".bg"
    else:
        return ".cg"

def bg_to_elem_string(cg):
    return fgb.BulgeGraph.to_dotbracket_string(cg)+"\n"+fgb.BulgeGraph.to_element_string(cg, with_numbers=True)

class OutFiletype:
    def __init__(self, write_fun, extention_fun, rna_type):
        self.convert = write_fun
        self.get_extention = extention_fun
        self.rna_type = rna_type

FILETYPES = {
    "forgi":OutFiletype(to_bg_or_cg_string, cg_or_bg_extention, "any"),
    "bpseq":OutFiletype(fgb.BulgeGraph.to_bpseq_string,lambda x: ".bpseq","cg"),
    "fasta":OutFiletype(fgb.BulgeGraph.to_fasta_string,lambda x: ".fa", "cg"),
    "dotbracket":OutFiletype(fgb.BulgeGraph.to_dotbracket_string, lambda x: ".dotbracket", "any"),
    "neato":OutFiletype(fgb.BulgeGraph.to_neato_string, lambda x: ".neato", "any"),
    "element_string":OutFiletype(bg_to_elem_string, lambda x: ".element_string", "any"),
}

parser = fuc.get_rna_input_parser("Convert RNA files between different file formats.", "+")
parser.add_argument("-T", "--target-type", default="forgi", choices=FILETYPES.keys(), type=str, help="The target file-type to convert into.")
parser.add_argument("--to-file", action="store_true",
                    help="Store the converted RNA in files instead of printing them to stdout. "
                         "The file-name will be the RNA's name (if present), otherwise 'rna001' etc.")
parser.add_argument("--filename", type=str,
                    help="If this is present, --to-file will automatically be true."
                         "A target filename (or path) without extention. "
                         "If it is a filename, use the given filename instead of the RNA's name. "
                         "If more than one input-RNA is present, appends automatically a increasing number."
                         "If it is a directory, create files in this directory.")
parser.add_argument("-f", "--force", action="store_true",
                    help="Overwrite files, if they already exist. Note: In case of race conditions, "
                         "files could be overwritten even if this flag is not provided.")



if __name__=="__main__":
    args = parser.parse_args()
    cgs = fuc.cgs_from_args(args, "+", rna_type=FILETYPES[args.target_type].rna_type)
    if args.filename:
        args.to_file = True
        if os.path.isdir(args.filename):
            directory = args.filename
            filename = None
            args.filename = None
        else:
            directory, filename = os.path.split(args.filename)
    else:
        filename = None
        directory = ""
    for i, cg in enumerate(cgs):
        if not args.to_file and i>0:
            print("\n\n========================================\n\n")
        target_str = FILETYPES[args.target_type].convert(cg)
        if args.to_file:
            if filename:
                if len(cgs)>0:
                    fn = filename+"{:03d}".format(i+1)
                else:
                    fn = filename
            else:
                if cg.name and cg.name != "untitled":
                    fn = cg.name
                else:
                    fn = "rna{:03d}".format(i+1)
            fn+=FILETYPES[args.target_type].get_extention(cg)
            fn = os.path.join(directory, fn)
        else:
            fn="-"
        with fuc.open_for_out(fn, args.force) as outfile:
            print(target_str, file=outfile)
            if fn !="-":
                log.info("file {} written".format(fn))
