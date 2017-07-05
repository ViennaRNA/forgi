from __future__ import print_function

import logging

import numpy as np

import forgi.utilities.commandline_utils as fuc
import forgi.graph.bulge_graph as fgb

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

CONVERSIONS = {
    "forgi":to_bg_or_cg_string,
    "bpseq":fgb.BulgeGraph.to_bpseq_string,
    "fasta":fgb.BulgeGraph.to_fasta_string,
    "dotbracket":fgb.BulgeGraph.to_dotbracket_string
}
EXTENTIONS = {
    "forgi": cg_or_bg_extention,
    "bpseq": lambda x: ".bpseq",
    "fasta": lambda x: ".fa",
    "dotbracket": lambda x: ".dotbracket"
}

parser = fuc.get_parser_any_cgs("Convert RNA files between different file formats.", "+")
parser.add_argument("-T", "--target-type", default="forgi", choices=CONVERSIONS.keys(), type=str, help="The target file-type to convert into.")
parser.add_argument("--to-file", action="store_true",
                    help="Store the converted RNA in files instead of printing them to stdout. "
                         "The file-name will be the RNA's name (if present), otherwise 'rna001' etc.")
parser.add_argument("--filename", type=str,
                    help="A target filename without extention. "
                         "Instead of using the RNA's name from the input, use the given filename. "
                         "If more than one input-RNA is poresent, appends automatically a increasing number.")
parser.add_argument("-f", "--force", action="store_true",
                    help="Overwrite files, if they already exist. Note: In case of race conditions, "
                         "files could be overwritten even if this flag is not provided.")



if __name__=="__main__":
    args = parser.parse_args()
    cgs = fuc.parse_any_cgs(args, "+")
    for i, cg in enumerate(cgs):
        if not args.to_file and i>0:
            print("\n\n========================================\n\n")
        target_str = CONVERSIONS[args.target_type](cg)
        if args.to_file:
            if args.filename:
                if len(cgs)>0:
                    filename = args.filename+"{:03d}".format(i+1)
                else:
                    filename = args.filename
            else:
                if cg.name:
                    filename = cg.name
                else:
                    filename = "rna{:03d}".format(i+1)
            filename+=EXTENTIONS[args.target_type](cg)
        else:
            filename="-"
        with fuc.open_for_out(filename, args.force) as outfile:
            print(target_str, file=outfile)
            if filename !="-":
                log.info("file {} written".format(filename))
