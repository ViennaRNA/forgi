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

def get_parser_any_cgs(helptext, nargs = 1, rna_type = "any", enable_logging=True, parser_kwargs={}):
    parser = argparse.ArgumentParser(helptext, **parser_kwargs)
    if nargs==1:
        helptext = "One file containing an RNA.\n"
    elif isinstance(nargs, int) and nargs>1:
        helptext="{:d} files containing one RNA each.\n".format(nargs)
    elif nargs == "+":
        helptext="One or more files containing one or more RNAs each.\n"
    else:
        raise ValueError("get_parser_any_cgsg does not support nargs={}".format(nargs) )
    helptext+=" Supported file-formats are: \n  pdb files, forgi cg files,\n"
    if rna_type!="3d":
        helptext+=" forgi bg files, fasta and dotbracket files.\n "
    if rna_type=="any":
        helptext+=("Alternatively you can supply a dotbracket-string\n "
                   "(containing only the characters '.()[]{}&') from the commandline.\n")
    parser.add_argument("rna", nargs = nargs, type=str, help=helptext)
    parser.add_argument("--pseudoknots", action="store_true", help="Allow pseudoknots when extracting the structure from PDB files.")
    parser.add_argument("--chain", type=str,
                        help="When reading pdb-files: Only extract the given chain.")
    parser.add_argument("--pdb-secondary-structure", type=str, default="",
                        help="When reading a single chain from a pdb-files: \nEnforce the given secondary structure (as dotbracket string).\n (This only works, if --chain is given!)")

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

def load_rna(filename, rna_type="any", allow_many=True, pdb_chain=None, pbd_remove_pk=True, pdb_dotbracket=""):
    """
    :param rna_type: One of "any", "cg" and "3d"
                     "any": Return either BulgeGraph or CoarseGrainRNA objekte, depending on the input format
                     "cg": Always convert to CoarseGrainRNA objects, even if they have no 3D information
                     "3d": Return CoarseGrainRNA objects, if the file contains 3D information, raise an error otherwise
    :param allow_many: If True, return a list. If False raise an error, if more than one RNA is present.
    :param pdb_chain: Extract the given chain from the file. Only applicable if filename corresponds to a pdb file
    :param pdb_remove_pk: Detect pseudoknot-free structures from the pdb.
    :param odb_dotbracket: Only applicable, if filename corresponds to a pdb file and pdb_chain is given.
    """
    with open(filename) as rnafile:
        filetype = sniff_filetype(rnafile)
    if filetype=="forgi":
        cg = ftmc.CoarseGrainRNA(filename)
        if rna_type=="3d" and np.isnan(cg.coords).any():
            raise ValueError("File {} does not contain all 3D coordinates!".format(filename))
        if allow_many:
            return [cg]
        else:
            return cg
    elif filetype=="pdb":
        if pdb_chain:
            cgs = [ftmc.load_cg_from_pdb(filename, chain_id=pdb_chain,
                                                 remove_pseudoknots=pbd_remove_pk and not pdb_dotbracket,
                                                 secondary_structure=pdb_dotbracket)]
        else:
            cgs = ftmc.connected_cgs_from_pdb(filename, remove_pseudoknots = not args.pseudoknots)
        if allow_many:
            return cgs
        else:
            if len(cgs)>1:
                raise ValueError("More than one connected RNA component in pdb file {}.".format(filename))
            return cgs[0]
    elif filetype=="mmcif":
        raise ValueError("MMCIF files are not yet supported.")
    elif filetype=="bpseq":
        if rna_type=="3d":
            raise ValueError("bpseq file {} is not supported. We need 3D coordinates!".format(filename))
        bg = fgb.BulgeGraph()
        with open(filename, 'r') as f:
            text = f.read()
            try:
                int(text[0])
            except ValueError:
                i=text.find("\n1 ")
                text=text[i+1:]
        bg.from_bpseq_str(text)
        if rna_type=="cg":
            bg = ftmc.from_bulge_graph(bg)
        if allow_many:
            return [cg]
        else:
            return cg
    elif filetype =="fasta" or filetype=="other":
        if rna_type=="3d":
            raise ValueError("Fasta(like) file {} is not supported. We need 3D coordinates!".format(filename))
        try:
            bgs = fgb.from_fasta(filename)
        except Exception as e:
            with log_to_exception(log, e):
                log.critical("Could not parse file %r.", filename)
                if filetype=="other":
                    log.critical("We assumed file %r to be some fasta-variant or dotbracket file, but an error occurred during parsing.", filename)
            raise
        if isinstance(bgs, fgb.BulgeGraph):
            bgs = [bgs]
        if rna_type=="cg":
            bgs = map(ftmc.from_bulge_graph, bgs)
        if allow_many:
            return bgs
        else:
            if len(bgs)>1:
                raise ValueError("More than one RNA found in fasta/ dotbracket file {}.".format(filename))
            return bgs[0]



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


def parse_any_cgs(args, nargs = 1, rna_type="cg", enable_logging=True):
    if enable_logging:
        logging.basicConfig()
        logging_exceptions.config_from_args(args)

    cg_rnas = []
    for rna in args.rna:
        if rna_type=="any" and all( c in ".()[]{}&" for c in rna):
            # A dotbracket-string was provided via the commandline
            log.info("Assuming RNA %s is a dotbracketstring and not a file.", rna)
            bg = fgb.from_fasta_text(rna)
            cg_rnas.append(bg)
        elif isinstance(nargs, int):
            cg_rnas.append(load_rna(rna, rna_type=rna_type, allow_many=False, pdb_chain=args.chain,
                                    pbd_remove_pk=not args.pseudoknots, pdb_dotbracket=args.pdb_secondary_structure))
        else:
            cg_rnas.extend(load_rna(rna, rna_type=rna_type, allow_many=True, pdb_chain=args.chain,
                                    pbd_remove_pk=not args.pseudoknots, pdb_dotbracket=pdb-secondary-structure))
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
