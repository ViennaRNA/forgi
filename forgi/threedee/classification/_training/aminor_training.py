from __future__ import print_function, absolute_import, division

import os
import argparse
try:
    from collections.abc import Set
except ImportError:
    from collections import Set
from collections import namedtuple, defaultdict
import warnings
import logging
import time
import subprocess

import numpy as np
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.metrics import confusion_matrix, f1_score

import pandas as pd

from logging_exceptions import log_to_exception
from .. import aminor as ftca
import forgi.utilities.stuff as fus
from ._parse_utils import ChainIdMappingParser
import forgi.graph.residue as fgr
import forgi.graph.bulge_graph as fgb


try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

log = logging.getLogger(__name__)


PDBIDS_LRSU = ["1HC8", "1MMS", "1MZP", "1VQO", "3jbu", "4IOA", "4LGT", "4QVI", "4v9p",
               "4w2g", "4wt8", "5D8H", "5imq", "5MMI", "5O60", "5V7Q", "3J79", "3J7Q", "5T2C"]
PDBIDS_SRSU = ["1G1X", "1I6U", "1KUQ", "2VQE", "4V19", "4v9o", "5lzd", "5mmm",
               "5o61", "5OOL", "5v93", "3J7A", "3J7P", "4UG0", "3JAM", "4P8Z",
               "4UG0", "4V5O", "4V88", "5FLX", "5LZS", "5OPT", "5XXU", "5XYI", "6AZ1"]


def create_geometry_file(outfile, cgs, cg_filenames, fr3d_result,
                         chain_id_mapping_dir, fr3d_query=""):
    ame, non_ame = from_fr3d_to_orientation(
        cgs, fr3d_result, chain_id_mapping_dir)
    # First of all, print the ame_geometries.csv file!
    with open(outfile, "w") as file_:
        print("# Generated on " + time.strftime("%d %b %Y %H:%M:%S %Z"), file=file_)
        print("# Generated from the following files: ", file=file_)
        for fn in cg_filenames:
            print("#   {}".format(fn), file=file_)
        print("# Version: {}".format(fus.get_version_string().strip()), file=file_)
        print("# fr3d_result = {}".format(fr3d_result), file=file_)
        if fr3d_query:
            print("# fr3d_query:", file=file_)
            for line in fr3d_query.splitlines():
                line = line.strip()
                print("#    " + line, file=file_)
        print("# cutoff_dist = {} A".format(ftca.CUTOFFDIST), file=file_)
        # HEADER
        print ("pdb_id loop_type dist angle1 angle2 is_interaction "
               "loop_sequence interaction score annotation loop_flexibility", file=file_)
        for entry in ame:
            print("{pdb_id} {loop_type} {dist} {angle1} "
                  "{angle2} {is_interaction} {loop_sequence} "
                  "{loop_name}-{stem_name} {score} "
                  "\"{annotation}\" {loop_flexibility}".format(is_interaction=True,
                                                               **entry._asdict()), file=file_)
        for entry in non_ame:
            print("{pdb_id} {loop_type} {dist} {angle1} "
                  "{angle2} {is_interaction} {loop_sequence} "
                  "{loop_name}-{stem_name} {score} "
                  "\"{annotation}\" {loop_flexibility}".format(is_interaction=False,
                                                               **entry._asdict()), file=file_)


def tune_model(geometry_file, trainset, testset):
    if not trainset:
        trainset = PDBIDS_LRSU
    if not testset:
        testset = PDBIDS_SRSU
    df = pd.read_csv(geometry_file, comment="#", sep=" ")
    all_params = {}
    for loop_type in "hi":
        pI = calculate_pI(df, loop_type)
        # sum(labels)/(sum(labels)+len(df[mask_non_fred]))
        log.info("PI %s", pI)
        all_params[loop_type] = find_hyperparameters(loop_type, df,
                                                     pI, trainset, testset)
    return all_params

################################################################################
# Validate Model
################################################################################


def calculate_pI(df, loop_type):
    mask_ame, mask_non_ame, mask_non_fred = ftca._get_masks(df, loop_type)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df["loop_name"] = df.pdb_id + \
            df.interaction.str.split("-").apply(lambda x: x[0])
    negative = len(df[mask_non_fred].drop_duplicates(["loop_name"]))
    positive = len(df[mask_ame].drop_duplicates("loop_name"))
    return positive / (negative + positive)


def mask_data(df, pdbids):
    mask = df.pdb_id.str.upper().str.contains("testestest")
    for pdb_id in pdbids:
        mask |= df.pdb_id.str.upper().str.contains(pdb_id.upper())
    return df[mask]


def find_hyperparameters(loop_type, df, pI, trainset, testset):
    log.info("Now searching for bandwidth")
    df_train = mask_data(df, trainset)
    df_test = mask_data(df, testset)
    if len(df_train) == 0 or len(df_test) == 0:
        raise ValueError("Either trainings- or testset are empty.")
    X_train, y_train = ftca.df_to_data_labels(df_train, loop_type)
    X_test, y_test = ftca.df_to_data_labels(df_test, loop_type)
    log.info("Length trainingsset: %s, length testset: %s",
             len(X_train), len(X_test))

    # We use the gaussian kernel to avoid certain artefacts of the linear
    # one and symmetric due to geometric considerations.
    clf = ftca.AMinorClassifier(p_I=pI, kernel="gaussian", symmetric=True)
    clf.fit(X_train, y_train)
    scores = []
    bandwidths = np.linspace(0.01, 1.0, 25)
    for bandwidth in bandwidths:
        log.debug("bandwidth %s", bandwidth)
        clf.set_params(bandwidth=bandwidth)
        score = clf.score(X_test, y_test)
        log.debug("score is %s", score)
        scores.append(score)
    best_i = scores.index(max(scores))
    best_bandwidth = bandwidths[best_i]
    clf.set_params(bandwidth=best_bandwidth)
    y_pred = clf.predict(X_test)
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    specificity = tn / (tn + fp)
    sensitivity = tp / (tp + fn)
    log.info("Classification report for loop type %s (using CV, bw=%s)\n:"
             "tp: %s, fp: %s, tn: %s, fn: %s, sens: %s, spec; %s",
             loop_type, best_bandwidth, tp, fp, tn, fn, sensitivity, specificity)
    if sensitivity < 0.8 or specificity < 0.8:
        warnings.warn("The crossvalidation score for loop tpype {} is "
                      "worse than expected: sensitivity: {},"
                      " specificity: {}.".format(loop_type, sensitivity, specificity))
    if sensitivity < 0.5 or specificity < 0.7:
        raise RuntimeError("Something went wrong with the training of the "
                           "classifier. The crossvalidation score is "
                           "terrible: sensitivity: {},"
                           " specificity: {}.".format(sensitivity, specificity))
    return clf.get_params()


################################################################################
# From FRED to TRAININGSDATA
################################################################################
_AMGeometry = namedtuple("AMGeometry",
                         ["pdb_id", "loop_name", "stem_name", "dist", "angle1",
                          "angle2", "loop_sequence", "score", "annotation",
                          "loop_flexibility"])


class AMGeometry(_AMGeometry):
    def _asdict(self):
        d = super(AMGeometry, self)._asdict()
        d["loop_type"] = self.loop_type
        return d

    @property
    def loop_type(self):
        return self.loop_name[0]


class AmeGeometrySet(Set):
    """
    Like a set, but use only pdb_id, loop_name and stem_name for
    equality, not the hash.

    If multiple geometries with the same key are added, use the lowest-scoring.
    """

    def __init__(self):
        self.geometries = {}

    @staticmethod
    def _geomerty_to_key(geo):
        return (geo.pdb_id, geo.loop_name, geo.stem_name)

    def add(self, geometry):
        key = self._geomerty_to_key(geometry)
        if key in self.geometries:
            # Smaller score is better
            if geometry.score > self.geometries[key].score:
                log.info("Duplicate geometry: %s has worse score "
                         "than %s", geometry, self.geometries[key])
                return
            else:
                log.info("Duplicate geometry: "
                         "%s replacing %s", geometry, self.geometries[key])
        self.geometries[key] = geometry

    def __contains__(self, geo):
        key = self._geomerty_to_key(geo)
        return key in self.geometries

    def __iter__(self):
        for g in self.geometries.values():
            yield g

    def __len__(self):
        return len(self.geometries)


def _chains_from_fr3d_field(value, pdb_id, mapping_directory):
    if len(value) == 3:
        chains = value
    elif len(value) == 6:
        chains = [value[:2], value[2:4], value[4:]]
        # Asterix as a special character for pdbs with length 1-chain-ids in bundles
        chains = list(map(lambda x: x.replace('*', ''), chains))
    else:
        raise ValueError(
            "PDB-id {}: Expecting either 3 or 6 letters for chain, found '{}'.".format(pdb_id, value))
    try:
        chain_id_mapping_file = os.path.join(
            mapping_directory, pdb_id.lower() + "-chain-id-mapping.txt")
        mapping = ChainIdMappingParser().parse(chain_id_mapping_file).mmcif2bundle
    except FileNotFoundError:
        return chains
    else:
        return [mapping[c] for c in chains]


def _parse_fred_line(line, all_cgs, current_annotation, mapping_directory):
    parts = line.split()
    if len(parts) < 10:
        return
    if parts[0] == "Filename":
        return
    pdb_id = parts[0]

    if parts[0] in all_cgs:
        chains = _chains_from_fr3d_field(parts[8], parts[0], mapping_directory)
        if isinstance(chains, list):  # A pdb bundle:
            if not all(c[0] == chains[0][0] for c in chains):
                log.info("Found an Interaction between multiple "
                         "files in pdb-bundle: {} {}. "
                         "Ignoring it!".format(parts[0], [c[0] for c in chains]))
                return
            pdb_name = chains[0][0].split(".")[0]
            chains = [c[1] for c in chains]
            log.debug(
                "For bundle: pdb_name = %s, converted-chains = %s", pdb_name, chains)
        else:
            pdb_name = parts[0]
        cgs = all_cgs[parts[0]]
        if len(cgs) == 1:
            cg, = cgs
        else:
            # First, filter by pdb-bundle file.
            cgs = [cg for cg in cgs if pdb_name in cg.name]
            log.debug("cgs for pbd-name %s are %s",
                      pdb_name, [c.name for c in cgs])
            # Then for multiple chains, search based on resid.
            a_res = _safe_resid_from_chain_res(
                chain=chains[0], residue=parts[3])
            if a_res is None:
                warnings.warn("Cannot create resid for {}".format(chains[0]))
                return
            for cg in cgs:
                # Find the first matching cg.
                if a_res in cg.seq._seqids:
                    break
            else:
                warnings.warn(
                    "No CoarseGrainRNA found among those with PDB-ID {} that contains resid {}".format(parts[0], a_res))
                return
    else:
        warnings.warn(
            "No CoarseGrainRNA found for FR3D annotation with PDB-ID {}. Possible PDB-IDs are {}".format(parts[0], all_cgs.keys()))
        return
    log.debug("FR3D line refers to cg %s", cg.name)
    # We now have the CoarseGrainRNA. Get the interacting cg-elements.
    nums = []
    for i in range(3):
        seq_id = _safe_resid_from_chain_res(
            chain=chains[i], residue=parts[3 + 2 * i])
        if seq_id is None:
            warnings.warn("Cannot create resid for {}".format(chains[i]))
            return
        try:
            nums.append(cg.seq_id_to_pos(seq_id))
        except ValueError:
            if seq_id.chain not in cg.chains:
                warnings.warn("Chain {!r} is not part of the cg {}. Required chains: {}. Available "
                              "cgs: {}".format(seq_id.chain, cg.name, chains,
                                               list(map(lambda x: x.name, all_cgs[parts[0]]))))
                return
            else:
                log.error("%s %s", cg.seq, type(cg.seq))
                log.error("All seq_ids=%s", cg.seq._seqids)
                raise
    nodes = list(map(lambda x: cg.get_node_from_residue_num(x), nums))
    if nodes[1] != nodes[2] or nodes[1][0] != "s":
        log.info("Parse_fred: No canonical stem: {} != {}.".format(
            nodes[1], nodes[2]))
        return  # only look at canonical stems.
    if nodes[0][0] == "s":
        log.info(
            "Parse_fred: Stem-Stem A-Minor not allowed: {} -> {}.".format(nodes[0], nodes[1]))
        return  # No A-Minor between two stems.
    if nodes[0] in cg.edges[nodes[1]]:
        log.info(
            "Parse_fred:  A-Minor between adjacent elements not allowed: {} -> {}.".format(nodes[0], nodes[1]))
        return  # Only look at not connected elements
    dist, angle1, angle2 = ftca.get_relative_orientation(
        cg, nodes[0], nodes[1])
    if np.isnan(angle1 + angle2 + dist):
        warnings.warn(
            "Cannot get relative orientation. Zero-length element {}".format(nodes[0]))
        return
    if nodes[0][0] == "i":
        flexibility = ftca.get_loop_flexibility(cg, nodes[0])
    else:
        flexibility = 1
    return (AMGeometry(cg.name, nodes[0], nodes[1], dist, angle1, angle2,
                       "&".join(cg.get_define_seq_str(nodes[0])),
                       float(parts[1]), current_annotation,
                       flexibility))


def _safe_resid_from_chain_res(chain, residue):
    try:
        return fgr.resid_from_str(str("{}:{}".format(chain, residue)))
    except ValueError as e:
        if residue.isdigit():
            with log_to_exception(log, e):
                log.error("Chain is '{}', res is '{}'".format(chain, residue))
            raise
        else:
            warnings.warn("Illegal residue number: '{}'.".format(residue))
            return


def parse_fred(cutoff_dist, all_cgs, fr3d_out, chain_id_mapping_dir):
    """
    Used by generate_target_distribution.

    :param all_cgs: A dictionary PDB_id (4-letter): [CoarseGrainRNA, ...]
    :param fr3d_out: A file opened for reading.

    :returns: A set of AMinor Geometries and a list of PDB-IDs, for which at
              least one line did not lead to an Aminor Geometry
    """
    with warnings.catch_warnings():
        warnings.simplefilter('always', UserWarning, append=True)
        geometries = AmeGeometrySet()
        skipped = 0
        #: What type of AMinor interactions.
        #: Comments in the form "# AMinor 0" have to be added manually
        #: when copying the FR3D-output to a file.They should preced the
        #: actual FR3D-output to which they apply.
        current_annotation = "?"
        for line in fr3d_out:
            line = line.strip()
            if not line:  # Empty line
                continue
            if line.startswith("# AMinor"):
                current_annotation = line[9:]
            elif line.startswith("#"):
                current_annotation = "?"
            log.debug("Line '%s'.read", line)
            geometry = _parse_fred_line(
                line, all_cgs, current_annotation, chain_id_mapping_dir)
            if geometry is None:
                if not (line.startswith("Filename") or line.startswith("#")):
                    skipped += 1
                    log.info("Skipping line {!r}".format(line))
            elif geometry.dist > cutoff_dist:
                log.debug("Skipping because of %f > %f (=cutoff dist): %r",
                          geometry.dist, cutoff_dist, line)
            elif "A" not in geometry.loop_sequence:
                warnings.warn("No adenine in loop %r for line %r",
                              geometry.loop_name, line)
            else:
                geometries.add(geometry)
    return geometries, skipped


def from_fr3d_to_orientation(cgs, fr3d_outfile, chain_id_mapping_dir):
    """
    Generate the trainings data for the AMinor classification
    from a list of cgs and the FR3D results of a query against
    the corresponding pdbs.

    ..note::

        The FR3D query needs to contain 3 nucleotides where the first is the Adenine and
        the second and third form a basepair in a stem.

    It needs FR3D output like this::

        Filename Discrepancy       1       2       3 Cha   1-2   1-3   2-3 Con   1-2   1-3   2-3   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-3
            1S72      0.0340  A  104  A  957  U 1009 900 ----  ----   cWW  AAA  1909  1885    24                                                                             -     -     -

    :param cgs: A list of CoarseGrainRNAs
                .. note::

                    the first 4 characters of the CoarseGrainRNA's name must match the
                    pdb-id as reported by FR3D.

    :param fr3d_out: A file containing the annotations found by FR3D.
    :param orientation_outfile: Where to store the list of relative orientations of
                                Adenines and the Stems in A-Minor interactions and
                                in the background. Relative to "fess/"
    :param fr3d_query: Describe the FR3D-query that you used in a string.
                       It will be added as comment to the output file.
    :param chain_id_mapping_dir: A dictionary containing the
                      chain-id-mapping files for pdb-bundles.
    """
    # Create lookup-dictionary from cgs.
    all_cgs = defaultdict(list)
    for cg in cgs:
        all_cgs[cg.name[:4].upper()].append(cg)
    # Read the FR3D output
    with open(fr3d_outfile) as f:
        aminor_geometries, num_skipped = parse_fred(
            ftca.CUTOFFDIST, all_cgs, f, chain_id_mapping_dir)
    log.info("%d entries skipped during FR3D parsing, %s entries retained",
             num_skipped, len(aminor_geometries))
    if len(aminor_geometries) == 0:
        raise ValueError(
            "No A-Minor geometries found. Is the FR3D output file correct?")
    non_ame_geometries = _enumerate_background_geometries(
        all_cgs, ftca.CUTOFFDIST, aminor_geometries)
    return aminor_geometries, non_ame_geometries


def _enumerate_background_geometries(all_cgs, cutoff_dist, aminor_geometries):
    """
    :param all_cgs: A dictionary {PDBID: [ cg1, cg2, ...]}
    """
    non_ame_geometries = set()
    for pdb_id, curr_cgs in all_cgs.items():
        for cg in curr_cgs:
            try:
                for loop in cg.defines:
                    if loop[0] == "s":
                        continue
                    if loop in cg.incomplete_elements or loop in cg.interacting_elements:
                        continue
                    for stem in cg.stem_iterator():
                        if loop in cg.edges[stem]:
                            continue
                        if stem in cg.incomplete_elements or stem in cg.interacting_elements:
                            continue
                        dist, angle1, angle2 = ftca.get_relative_orientation(
                            cg, loop, stem)
                        if loop[0] == "i":
                            flexibility = ftca.get_loop_flexibility(cg, loop)
                        else:
                            flexibility = 1
                        if not np.isnan(dist + angle1 + angle2) and dist <= cutoff_dist:
                            geometry = AMGeometry(cg.name, loop, stem, dist,
                                                  angle1, angle2,
                                                  "&".join(
                                                      cg.get_define_seq_str(loop)),
                                                  1000, "no_interaction",
                                                  flexibility)
                            if geometry in aminor_geometries:
                                log.info(
                                    "Geometry %s is in aminor_geometries", geometry)
                            else:
                                non_ame_geometries.add(geometry)
            except BaseException as e:
                with log_to_exception(log, e):
                    log.error(
                        "An Error occurred during processing of cg: %s", cg.name)
                raise
    log.error("%s non_ame geometries found", len(non_ame_geometries))
    return non_ame_geometries
