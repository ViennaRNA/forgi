from __future__ import print_function

import sys
import argparse
import logging
import os.path
import contextlib
import warnings
import textwrap
import copy

import numpy as np

from .exceptions import GraphConstructionError
from logging_exceptions import log_to_exception
import logging_exceptions


import forgi.graph.bulge_graph as fgb
import forgi.graph.sequence as fgs
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
        pdb_input_group.add_argument("--pdb-allow-www-query", action="store_true",
                                     help="Usually, if modified residues/ ligand with an unknown 3-letter code are encountered in PDB files, "
                                          "they are removed from the chain and a log message of severity INFO is issued. "
                                          "With this option, we first try to query the PDBeChem database for the 3-letter code, to "
                                          "see whether or not it is a modified residue that can be converted to its standard parent "
                                          " and should be part of the chain.")

    if enable_logging:
        verbosity_group = parser.add_argument_group(
            "Control verbosity of logging output")
        logging_exceptions.update_parser(verbosity_group)
    return parser


def cgs_from_args(args, rna_type="any", enable_logging=True,
                  return_filenames=False, skip_errors=False):
    """
    Given an Namespace from argparse, return a list of CoarseGrainRNA objects
    (or BulgeGraph objects).

    :param args: A argparse Namespace.
    :param rna_type: See documentation of load_rna
    :param enable_logging: In addition to loading the CoarseGrainRNA objects,
                    call logging.basicConfig() and use args.verbose
                    and args.debug
                    to set the level for different loggers.
    :param return_filenames: Instead of returning a list of rnas,
                    return a tuple `rnas, filenames`,
                    where filenames is a list of filenames of equal length as
                    the returned list of CoarseGrainRNA objects and with
                    corresponding order.
                    If a file contains two seperate RNA molecules, the
                    filename will thus be found twice in the list of filenames.
    :param skip_errors: Boolean. Log GraphConstructionErrors and continue with
                    the next filename instead of letting the error propagate.

    Usage::

        parser = get_rna_input_parser()
        args = parser.parse_args()
        rnas = cgs_from_args(args)

    """
    if enable_logging:
        logging.basicConfig(
            format="%(levelname)s:%(name)s.%(funcName)s[%(lineno)d]: %(message)s")
        logging_exceptions.config_from_args(args)
    cg_rnas = []
    filenames = []
    for rna in args.rna:
        log.debug("Load RNA %s", rna)
        if rna_type == "only_cg":
            args.chains = None
            args.pseudoknots = None
            args.pdb_secondary_structure = None
            args.pdb_annotation_tool = None
            args.pdb_allow_www_query = False
        if args.chains:
            load_chains = args.chains.split(",")
        else:
            load_chains = None
        try:
            cg_or_cgs = load_rna(rna, rna_type=rna_type, allow_many=True,
                                 pdb_chain=load_chains,
                                 pdb_remove_pk=not args.pseudoknots, pdb_dotbracket=args.pdb_secondary_structure,
                                 dissolve_length_one_stems=not args.keep_length_one_stems,
                                 pdb_annotation_tool=args.pdb_annotation_tool,
                                 pdb_allow_www_query=args.pdb_allow_www_query)
        except GraphConstructionError:
            if not skip_errors:
                log.error("An error occurred while loading the file %s", rna)
                raise
            else:
                log.exception("The PDB %s was skipped due to the following error", rna)
        else:
            cg_rnas.extend(cg_or_cgs)
            filenames.extend([rna] * len(cg_or_cgs))
    if return_filenames:
        return cg_rnas, filenames
    else:
        return cg_rnas


def sniff_filetype(file):
    line = next(file)
    # PDB
    if line.startswith("ATOM") or line.startswith("HEADER") or line.startswith("HETATM") or line.startswith("REMARK"):
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
             pdb_annotation_tool=None, pdb_allow_www_query=False):
    """
    :param rna_type: One of "any", and "3d" and "pdb"

                     *  "any": Return either BulgeGraph or CoarseGrainRNA object,
                               depending on the input format
                     *  "only_cg": Only accept cg-files.
                     *  "3d":  Return CoarseGrainRNA objects,
                               if the file contains 3D information,
                               raise an error otherwise
                     *   "pdb": only accept pdb files

    :param allow_many: If True, return a list. If False, return a single
                       CoarseGrainRNA object or raise a WrongFileFormat,
                       if more than one RNA is present.
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
        if rna_type in ["3d", "only_cg"] and not cg.coords.is_filled: # pylint: disable=E1101
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
                                               filetype=filetype, annotation_tool=pdb_annotation_tool,
                                               query_PDBeChem=pdb_allow_www_query)
        else:
            if pdb_dotbracket:
                raise ValueError(
                    "pdb_dotbracket requires a chain to be given to avoid ambiguity.")
            cgs = ftmc.CoarseGrainRNA.from_pdb(filename, remove_pseudoknots=pdb_remove_pk,
                                               dissolve_length_one_stems=dissolve_length_one_stems,
                                               filetype=filetype, annotation_tool=pdb_annotation_tool,
                                               query_PDBeChem=pdb_allow_www_query)
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
            for record in e.log:
                logging.getLogger(record.name).handle(record)
        print("Error of type {} occurred. Aborting.".format(type(e).__name__))
        print(e, file=sys.stderr)
        sys.exit(1)

def db_to_constraint(dbstri):
    """
    Convert a dotbracket string with dashes representing missing residues and
    dots representing unpaired nucleotides to a constraint that can be used
    with RNAfold.

    This means that unpaired nucleotides are forced unpaired and missing
    residues are free to form basepairs.
    """
    dbstri = dbstri.replace(".", "x").replace("-",".")
    return ''.join( "x" if c not in "()." else c for c in dbstri )

def insert_pk_into_stru(refolded, with_pk):
    """
    Replaces '.' characters in refolded with pkseudoknots from with_pk where possible.

    :param with_pk: A dotbracket string, possibly with '-' characters
    :param refolded: The same dbstring as with_pk, but without pseudoknots and with '-' replaced by one of '().
    """
    out=[]
    for i,c in enumerate(refolded):
        if c != with_pk[i]:
            if c==".":
                if with_pk[i]=="-":
                    out.append(c)
                else:
                    out.append(with_pk[i])
            else:
                assert with_pk[i]=="-"
                out.append(c)
        else:
            out.append(c)
    return "".join(out)

def with_missing_refolded(cg):
    """
    This requires the ViennaRNA package.

    Uses RNAfold
    """
    domains = cg.get_domains()
    fasta_lines = cg.to_fasta_string(include_missing=True).splitlines()
    assert len(fasta_lines)==3
    constraint = db_to_constraint(fasta_lines[2])
    try:
        import RNA
    except ImportError as e:
        raise ImportError(str(e)+"\nPlease install the ViennaRNA package to fold missing residues.")
    fc = RNA.fold_compound(fasta_lines[1], RNA.md())
    fc.hc_add_from_db(constraint)
    stru, energy = fc.mfe()

    stru=insert_pk_into_stru(stru, fasta_lines[2])
    log.info("Energy of constraint secondary structyure is %s", energy)
    new_fasta="\n".join([fasta_lines[0], fasta_lines[1], stru])
    new_cg, = ftmc.CoarseGrainRNA.from_fasta_text(new_fasta)
    new_cg._seq = fgs.missing_to_normal(cg.seq)
    # TODO: Copy coordinates, but take potentially changed 2D structure into account.

    #new_cg.sampled = copy.deepcopy(cg.sampled)
    #new_cg.longrange = copy.deepcopy(cg.longrange)
    #new_cg.interacting_residues = copy.deepcopy(cg.interacting_residues)
    #new_cg.vposs = {k:copy.deepcopy(v) for k,v in cg.vposs.items() if k[0]!="s"}

    return new_cg
