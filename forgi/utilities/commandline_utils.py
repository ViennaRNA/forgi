import sys
import argparse
import logging
import os.path
import contextlib

import numpy as np

from logging_exceptions import log_to_exception
import logging_exceptions


import forgi.graph.bulge_graph as fgb
import forgi.threedee.model.coarse_grain as ftmc

log = logging.getLogger(__name__)

def get_parser_any_cgs(helptext, nargs = 1, allow_2d = True, enable_logging=True):
    parser = argparse.ArgumentParser(helptext)
    if nargs==1:
        helptext = "One file containing an RNA."
    elif isinstance(nargs, int) and nargs>1:
        helptext="{:d} files containing one RNA each.".format(nargs)
    elif nargs == "+":
        helptext="One or more files containing one or more RNAs each."
    else:
        raise ValueError("get_parser_any_cgsg does not support nargs={}".format(nargs) )
    helptext+=" Supported file-formats are: pdb files, forgi cg files"
    if allow_2d:
        helptext+=(" forgi bg files, fasta and dotbracket files. "
                  "Alternatively you can supply a dotbracket-string "
                  "(containing only the characters '.()[]{}&') from the commandline ")
    parser.add_argument("rna", nargs = nargs, type=str, help=helptext)
    parser.add_argument("-p", "--pseudoknots", action="store_true", help="Allow pseudoknots when extracting the structure from PDB files.")
    parser.add_argument("-c", "--chain", type=str,
                        help="When reading pdb-files: Only extract the given chain.")
    parser.add_argument("--pdb-secondary-structure", type=str, default="",
                        help="When reading a single chain from a pdb-files: Enforce the given secondary structure (as dotbracket string). (This only works, if --chain is given!)")

    if enable_logging:
        logging_exceptions.update_parser(parser)
    return parser
def sniff_filetype(file):
    line = next(file)
    if line.startswith("ATOM") or line.startswith("HEADER"):
        return "pdb"
    line=line.strip()
    #We allow comments in all files except PDB files.
    while line.startswith("#"):
        line=next(file).strip()
    if line.startswith("name"):
        return "forgi"
    elif line.startswith("_") or line=="loop_":
        return "mmcif"
    elif line.startswith(">"):
        return "fasta"
    elif line[0].isdigit():
        return "bpseq"
    else:
        while True:
            line = next(file).strip()
            if line.startswith("#"):
                continue
            if line[0].isdigit():
                try:
                    d1, nt, d2 = line.split()
                    if d1.isdigit() and nt in "AUGCTIaugcti" and d2.isdigit():
                        return "bpseq"
                except (TypeError, ValueError):
                    pass
                break
        return "other"

def parse_any_cgs(args, nargs = 1, allow_2d=True, enable_logging=True):
    if enable_logging:
        logging.basicConfig()
        logging_exceptions.config_from_args(args)

    cg_rnas = []
    for rna in args.rna:
        if allow_2d and all( c in ".()[]{}&" for c in rna):
            # A dotbracket-string was provided via the commandline
            log.info("Assuming RNA %s is a dotbracketstring and not a file.", rna)
            cg_rnas.append(fgb.from_dotbracket(rna))
        else:
            # We do not look at the file extention, but instead at the file's content.
            with open(rna) as f:
                mode = sniff_filetype(f)
            if mode=="forgi":
                cg = fgb.BulgeGraph(rna)
                if not allow_2d and np.isnan(cg.coords).any():
                    raise ValueError("File {} does not contain all 3D coordinates!".format(rna))
                cg_rnas.append(cg)
            elif mode=="pdb":
                if args.chain:
                    cgs = [ftmc.load_cg_from_pdb(rna, chain_id=args.chain,
                                                 remove_pseudoknots=(not args.pseudoknots) and (not args.pdb_secondary_structure),
                                                 secondary_structure=args.pdb_secondary_structure)]
                else:
                    cgs = ftmc.connected_cgs_from_pdb(rna, remove_pseudoknots = not args.pseudoknots)
                if isinstance(nargs, int) and len(cgs)>1: #If nargs='+', we allow multiple RNAs per file.
                    raise ValueError("Only 1 RNA is allowed per file. More than 1 connected components found in pdb_file {}".format(rna))
                cg_rnas += cgs
            elif mode=="mmcif":
                raise ValueError("MMCIF files are not yet supported.")
            elif mode=="bpseq":
                if not allow_2d:
                    raise ValueError("bpseq file {} is not supported. We need 3D coordinates!".format(rna))
                bg = fgb.BulgeGraph()
                with open(rna, 'r') as f:
                    text = f.read()
                    try:
                        int(text[0])
                    except ValueError:
                        i=text.find("\n1 ")
                        text=text[i+1:]
                bg.from_bpseq_str(text)
                cg_rnas.append(bg)
            elif mode =="fasta" or mode=="other":
                if not allow_2d:
                    raise ValueError("Fasta(like) file {} is not supported. We need 3D coordinates!".format(rna))
                with open(rna) as f:
                    fasta_str = f.read()
                    try:
                        bgs = fgb.from_fasta_text(fasta_str)
                    except Exception as e:
                        with log_to_exception(log, e):
                            log.critical("Could not parse file %r.", rna)
                            if mode=="other":
                                log.critical("We assumed file %r to be some fasta-variant or dotbracket file, but an error occurred during parsing.", rna)
                        raise
                if isinstance(bgs, fgb.BulgeGraph):
                    bgs = [bgs]
                if isinstance(nargs, int) and len(bgs)>1: #If nargs='+', we allow multiple RNAs per file.
                    raise ValueError("Only 1 RNA is allowed per file. More than 1 found in fasta file {}".format(rna))
                cg_rnas += bgs
    return cg_rnas


@contextlib.contextmanager
def open_for_out(filename=None, force=False):
    "From http://stackoverflow.com/a/17603000/5069869"
    if filename and filename != '-':
        if not force and os.path.isfile(filename):
            raise IOError("Cannot create file {}. File exists.".format(filename))
        fh = open(filename, 'w')
    else:
        fh = sys.stdout
    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()
