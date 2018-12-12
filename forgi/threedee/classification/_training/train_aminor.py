from __future__ import print_function, absolute_import, division

import os
import argparse
import warnings
import logging
import json

import forgi.utilities.commandline_utils as fuc
import forgi.threedee.classification._training.aminor_training as ftcta
log = logging.getLogger(__name__)


PDBIDS_LRSU = ["1HC8", "1MMS", "1MZP", "1VQO", "3jbu", "4IOA", "4LGT", "4QVI", "4v9p",
               "4w2g", "4wt8", "5D8H", "5imq", "5MMI", "5O60", "5V7Q", "3J79", "3J7Q", "5T2C"]
PDBIDS_SRSU = ["1G1X", "1I6U", "1KUQ", "2VQE", "4V19", "4v9o", "5lzd", "5mmm",
               "5o61", "5OOL", "5v93", "3J7A", "3J7P", "4UG0", "3JAM", "4P8Z",
               "4UG0", "4V5O", "4V88", "5FLX", "5LZS", "5OPT", "5XXU", "5XYI", "6AZ1"]

################################################################################
# Lets use this as a script
# Run as python -m forgi.threedee.classification._training.aminor_training
################################################################################


def main():
    parser = fuc.get_rna_input_parser("Train a classifier for A-Minior interactions.",
                                      nargs="+", rna_type="only_cg")
    parser.add_argument("--fr3d-result", type=str, required=True,
                        help="A file containing the FR3D output")
    parser.add_argument("--fr3d-query", type=str,
                        help="Add this string describing the FR3D query "
                             "as a comment to the trainingsdata-out file")
    parser.add_argument("--chain-id-mapping-dir", type=str,
                        help="If you use PDB-bundles, this directory "
                             "needs to hold all chain-id-mapping.txt files.")
    parser.add_argument("--trainingsdata-out", type=str,
                        default="forgi/threedee/data/aminor_geometries.csv",
                        help="File that will be written for the geometries "
                             "of interactions and non-interactions.")
    parser.add_argument("--model-params-out", type=str,
                        default="forgi/threedee/data/aminor_params.json",
                        help="File that will be written for the "
                             "model's hyper-parameters.")
    parser.add_argument("--test-set", type=str, help="':'-separated PDB-ids"
                                                     " for the test-set.")
    parser.add_argument("--train-set", type=str,
                        help="':'-separated PDB-ids for the train-set."
                             "Note: This is only used for cross-validation."
                             "The final model will be trained on all the data.")

    args = parser.parse_args()
    cgs, cg_filenames = fuc.cgs_from_args(args, rna_type="only_cg",
                                          return_filenames=True)
    ftcta.create_geometry_file(args.trainingsdata_out, cgs, cg_filenames,
                               args.fr3d_result, args.chain_id_mapping_dir,
                               args.fr3d_query)
    hyper_params = ftcta.tune_model(
        args.trainingsdata_out, args.train_set, args.test_set)
    with open(args.model_params_out, "w") as f:
        json.dump(hyper_params, f)


if __name__ == "__main__":
    main()
