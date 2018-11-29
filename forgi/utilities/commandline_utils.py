from __future__ import print_function

import sys
import argparse
import logging
import os.path
import contextlib
import warnings
import textwrap

import numpy as np

from .exceptions import GraphConstructionError
from logging_exceptions import log_to_exception
import logging_exceptions


import forgi.graph.bulge_graph as fgb
import forgi.threedee.model.coarse_grain as ftmc

log = logging.getLogger(__name__)


class WrongFileFormat(ValueError):
    """
    Error raised if the input file has the wrong file format or the file format could not be detected.
    """


def get_rna_input_parser(helptext, nargs=1, rna_type="any", enable_logging=True, parser_kwargs={}):
    parser = argparse.ArgumentParser(description=helptext, **parser_kwargs)
    if nargs == 1:
        helptext = "One file containing an RNA.\n"
    elif isinstance(nargs, int) and nargs > 1:
        helptext = "{:d} files containing one RNA each.\n".format(nargs)
    elif nargs == "+":
        helptext = "One or more files containing one or more RNAs each.\n"
    else:
        raise ValueError(
            "get_parser_any_cgsg does not support nargs={}".format(nargs))

    if rna_type != "only_cg":
        fileformats = ["pdb files"]
    else:
        fileformats = []
    if rna_type != "pdb":
        fileformats.append("forgi cg files")
        if rna_type != "3d" and rna_type != "only_cg":
            fileformats.append("forgi bg files")
            fileformats.append("fasta files")
    if rna_type == "any":
        fileformats.append("dotbracketfiles")
    n = len(fileformats)
    helptext += "\n".join(textwrap.wrap(
        "Supported Filetypes are: {}".format(", ".join(fileformats)), 55))
    if rna_type == "any":
        helptext += ("Alternatively you can supply a dotbracket-string\n "
                     "(containing only the characters '.()[]{}&') from the commandline.\n")
    parser.add_argument("rna", nargs=nargs, type=str, help=helptext)
    parser.add_argument("--keep-length-one-stems", action="store_true",
                        help="For all input formats except forgi bg/cg files,\n"
                             "this controlls whether stems of length one are \n"
                             "dissolved to unpaired regions (default) \n"
                             "or kept (if this option is present).\n"
                             "In  the case of input in forgi-format,\n"
                             "the RNA from the file is not modified.")
    if rna_type != "only_cg":
        pdb_input_group = parser.add_argument_group("Options for loading of PDB files",
                                                    description="These options only take effect, "
                                                    "if the input RNA is in pdb file format.")
        pdb_input_group.add_argument("--pseudoknots", action="store_true",
                                     help="Allow pseudoknots when extracting the structure\nfrom PDB files.")
        pdb_input_group.add_argument("--chains", type=str,
                                     help="When reading pdb-files: Only extract the given chain(s). Comma-seperated")
        pdb_input_group.add_argument("--pdb-secondary-structure", type=str, default="",
                                     help="When reading a single chain from a pdb-files: \nEnforce the secondary structure given as dotbracket\n string. (This only works, if --chain is given!)")
        pdb_input_group.add_argument("--pdb-annotation-tool", type=str, default=None,
                                     help="What program to use for detecting basepairs in PDB/ MMCIF structures."
                                     " This commandline option overrides the value in the config file (if present)."
                                     "If this is not present and no config-file is given, we try to detect the installed programs.")

    if enable_logging:
        verbosity_group = parser.add_argument_group(
            "Control verbosity of logging output")
        logging_exceptions.update_parser(verbosity_group)
    return parser


def cgs_from_args(args, nargs=1, rna_type="any", enable_logging=True, return_filenames=False, skip_errors=False):
    if enable_logging:
        logging.basicConfig(
            format="%(levelname)s:%(name)s.%(funcName)s[%(lineno)d]: %(message)s")
        logging_exceptions.config_from_args(args)
    cg_rnas = []
    filenames = []
    for rna in args.rna:
        if isinstance(nargs, int):
            allow_many = False
        else:
            allow_many = True
        if rna_type == "only_cg":
            args.chains = None
            args.pseudoknots = None
            args.pdb_secondary_structure = None
            args.pdb_annotation_tool = None
        if args.chains:
            load_chains = args.chains.split(",")
        else:
            load_chains = None
        try:
            cg_or_cgs = load_rna(rna, rna_type=rna_type, allow_many=allow_many,
                                 pdb_chain=load_chains,
                                 pdb_remove_pk=not args.pseudoknots, pdb_dotbracket=args.pdb_secondary_structure,
                                 dissolve_length_one_stems=not args.keep_length_one_stems,
                                 pdb_annotation_tool=args.pdb_annotation_tool)
        except GraphConstructionError:
            if not skip_errors:
                raise
            else:
                log.exception("The following PDB was skipped")
        if allow_many:
            cg_rnas.extend(cg_or_cgs)
            filenames.extend([rna] * len(cg_or_cgs))
        else:
            cg_rnas.append(cg_or_cgs)
            filenames.append(rna)
    if return_filenames:
        return cg_rnas, filenames
    else:
        return cg_rnas


def sniff_filetype(file):
    line = next(file)
    # PDB
    if line.startswith("ATOM") or line.startswith("HEADER") or line.startswith("HETATM"):
        return "pdb"
    line = line.strip()
    # We allow comments and empty lines in all files except PDB files.
    while not line or line.startswith("#"):
        line = next(file).strip()
    if line.startswith("name"):
        return "forgi"
    elif line.startswith("_") or line == "loop_":
        return "cif"
    elif line.startswith(">"):
        return "fasta"
    elif line[0].isdigit():
        return "bpseq"
    else:
        # Is it a bp-seq or mmcif file with header?
        try:
            while True:
                line = next(file).strip()
                if line.startswith("#"):
                    continue
                if line.startswith("_"):
                    return "cif"
                if line[0].isdigit():
                    try:
                        d1, nt, d2 = line.split()
                        if d1.isdigit() and nt in "AUGCTIaugcti" and d2.isdigit():
                            return "bpseq"
                    except (TypeError, ValueError):
                        pass
                    break
        except StopIteration:
            pass
        return "other"


def load_rna(filename, rna_type="any", allow_many=True, pdb_chain=None,
             pdb_remove_pk=True, pdb_dotbracket="",
             dissolve_length_one_stems=True,
             pdb_annotation_tool=None):
    """
    :param rna_type: One of "any", and "3d" and "pdb"

                     *  "any": Return either BulgeGraph or CoarseGrainRNA objekte,
                               depending on the input format
                     *  "only_cg": Only accept cg-files.
                     *  "3d":  Return CoarseGrainRNA objects,
                               if the file contains 3D information,
                               raise an error otherwise
                     *   "pdb": only accept pdb files

    :param allow_many: If True, return a list. If False raise an error, if more than one RNA is present.
    :param pdb_chain: Extract the given chain from the file.
                      Only applicable if filename corresponds to a pdb file
    :param pdb_remove_pk: Detect pseudoknot-free structures from the pdb.
    :param pdb_dotbracket: Only applicable, if filename corresponds to a pdb file and pdb_chain is given.
    :param dissolve_length_one_stems: Ignored if input is in forgi bg/cg format.
    :param pdb_annotation_tool: Use DSSR, MC-Annotate or forgi heuristic for
                    basepair-detection in PDB/MMCIF files (None for auto-detect).
                    Ignored for other file-types.

    :retuns: A list of RNAs or a single RNA
    """
    # Is filename a dotbracket string and not a filename?
    if all(c in ".()[]{}&" for c in filename):
        # A dotbracket-string was provided via the commandline
        if not rna_type == "any":
            warnings.warn("Cannot treat '{}' as dotbracket string, since we need a sequence. "
                          "Trying to treat it as a filename instead...".format(filename))
        else:
            log.info(
                "Assuming RNA %s is a dotbracketstring and not a file.", filename)
            bg = fgb.BulgeGraph.from_dotbracket(
                filename, dissolve_length_one_stems=dissolve_length_one_stems)
            if allow_many:
                return [bg]
            else:
                return bg
    with open(filename) as rnafile:
        filetype = sniff_filetype(rnafile)
    if rna_type == "pdb" and filetype not in ["pdb", "cif"]:
        raise WrongFileFormat(
            "Only PDB files (*.pdb/.cif) are accepted, but file {} has type {}.".format(filename, filetype))
    if rna_type == "only_cg" and filetype != "forgi":
        raise WrongFileFormat(
            "Only forgi cg files are accepted, but file {} has type {}.".format(filename, filetype))
    if filetype == "forgi":
        cg = ftmc.CoarseGrainRNA.from_bg_file(filename)
        if rna_type in ["3d", "only_cg"] and not cg.coords.is_filled:
            raise WrongFileFormat(
                "File {} does not contain all 3D coordinates!".format(filename))
        if allow_many:
            return [cg]
        else:
            return cg
    elif filetype == "pdb" or filetype == "cif":
        if pdb_chain:
            cgs = ftmc.CoarseGrainRNA.from_pdb(filename, load_chains=pdb_chain,
                                               remove_pseudoknots=pdb_remove_pk and not pdb_dotbracket,
                                               secondary_structure=pdb_dotbracket,
                                               dissolve_length_one_stems=dissolve_length_one_stems,
                                               filetype=filetype, annotation_tool=pdb_annotation_tool)
        else:
            if pdb_dotbracket:
                raise ValueError(
                    "pdb_dotbracket requires a chain to be given to avoid ambiguity.")
            cgs = ftmc.CoarseGrainRNA.from_pdb(filename, remove_pseudoknots=pdb_remove_pk,
                                               dissolve_length_one_stems=dissolve_length_one_stems,
                                               filetype=filetype, annotation_tool=pdb_annotation_tool)
        if allow_many:
            return cgs
        else:
            if len(cgs) > 1:
                raise WrongFileFormat("More than one connected RNA component in pdb file {}: {}".format(
                    filename, [cg.name for cg in cgs]))
            return cgs[0]
    # elif filetype=="mmcif":
    #    raise WrongFileFormat("MMCIF files are not yet supported.")
    elif filetype == "bpseq":
        if rna_type == "3d":
            raise WrongFileFormat(
                "bpseq file {} is not supported. We need 3D coordinates!".format(filename))
        with open(filename, 'r') as f:
            text = f.read()
            try:
                int(text[0])
            except ValueError:
                i = text.find("\n1 ")
                text = text[i + 1:]
        bg = ftmc.CoarseGrainRNA.from_bpseq_str(
            text, dissolve_length_one_stems=dissolve_length_one_stems)
        if allow_many:
            return [bg]
        else:
            return bg
    elif filetype == "fasta" or filetype == "other":
        if rna_type == "3d":
            raise WrongFileFormat(
                "Fasta(like) file {} is not supported. We need 3D coordinates!".format(filename))
        try:
            bgs = ftmc.CoarseGrainRNA.from_fasta(
                filename, dissolve_length_one_stems=dissolve_length_one_stems)
        except Exception as e:
            with log_to_exception(log, e):
                log.critical("Could not parse file %r.", filename)
                if filetype == "other":
                    log.critical(
                        "We assumed file %r to be some fasta-variant or dotbracket file, but an error occurred during parsing.", filename)
            raise
        if allow_many:
            return bgs
        else:
            if len(bgs) > 1:
                raise WrongFileFormat(
                    "More than one RNA found in fasta/ dotbracket file {}.".format(filename))
            return bgs[0]


@contextlib.contextmanager
def open_for_out(filename=None, force=False):
    "From http://stackoverflow.com/a/17603000/5069869"
    if filename and filename != '-':
        if not force and os.path.isfile(filename):
            raise IOError(
                "Cannot create file {}. File exists.".format(filename))
        fh = open(filename, 'w')
    else:
        fh = sys.stdout
    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


@contextlib.contextmanager
def hide_traceback(error_class=WrongFileFormat):
    """
    A context manager that catches errors of the specified class and outputs
    only the error message before calling sys.exit(1).
    """
    try:
        yield
    except error_class as e:
        # Integration with logging_exceptions
        if hasattr(e, "log"):
            for record in exception.log:
                logging.getLogger(record.name).handle(record)
        print("Error of type {} occurred. Aborting.".format(type(e).__name__))
        print(e, file=sys.stderr)
        sys.exit(1)
