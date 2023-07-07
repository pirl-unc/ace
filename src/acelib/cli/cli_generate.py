# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
The purpose of this python3 script is to create parser
and run ACE 'generate' command.
"""


import math
import pandas as pd
import pkg_resources
import random
import torch
from ortools.sat.python import cp_model
from transformers import BertModel, BertTokenizer
from transformers import AutoTokenizer, AutoModelForMaskedLM
from ..constants import *
from ..default_parameters import *
from ..elispot import ELISpot
from ..logger import get_logger
from ..main import run_ace_sat_solver, run_ace_golfy
from ..utilities import convert_golfy_results, split_peptides
from ..sequence_features import AceNeuralEngine


logger = get_logger(__name__)


def add_ace_generate_arg_parser(sub_parsers):
    """
    Adds 'generate' parser.

    Parameters
    ----------
    sub_parsers     :   argparse.ArgumentParser subparsers object.

    Returns
    -------
    sub_parsers     :   argparse.ArgumentParser subparsers object.
    """
    parser = sub_parsers.add_parser(
        'generate',
        help='Generates an ELISpot experiment configuration.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required_mutually_exclusive = parser_required.add_mutually_exclusive_group(required=True)
    parser_required_mutually_exclusive.add_argument(
        "--num-peptides",
        dest="num_peptides",
        type=int,
        help="Total number of peptides."
    )
    parser_required_mutually_exclusive.add_argument(
        "--peptides-csv-file",
        dest="peptides_csv_file",
        type=str,
        help="Peptides CSV file with the following expected columns: "
             "'peptide_id', 'peptide_sequence'. Please note that only "
             "either this parameter or '--num-peptides' can be supplied."
    )
    parser_required.add_argument(
        "--num-peptides-per-pool",
        dest="num_peptides_per_pool",
        type=int,
        required=True,
        help="Number of peptides per pool (i.e. well)."
    )
    parser_required.add_argument(
        "--num-coverage",
        dest="num_coverage",
        type=int,
        required=True,
        help="Total coverage (i.e. number of peptide replicates)."
    )
    parser_required.add_argument(
        "--output-csv-file",
        dest="output_csv_file",
        type=str,
        required=True,
        help="Output CSV file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--mode",
        dest="mode",
        type=str,
        default=GenerateModes.GOLFY,
        choices=GenerateModes.ALL,
        required=False,
        help="Configuration generation mode (default: %s)." % GenerateModes.GOLFY
    )
    parser_optional.add_argument(
        "--assign-well-ids",
        dest="assign_well_ids",
        type=bool,
        default=True,
        required=False,
        help="If true, assigns plate and well IDs for each pool ID (default: true)."
    )
    parser_optional.add_argument(
        "--plate-type",
        dest="plate_type",
        type=str,
        default=PlateTypes.PLATE_96_WELLS,
        required=False,
        help="Plate type. Allowed values: %s (default: %s)." %
             (', '.join(PlateTypes.ALL), PlateTypes.PLATE_96_WELLS)
    )
    parser_optional.add_argument(
        "--trained-model",
        dest="trained_model",
        type=str,
        default=TrainedModels.MODEL_3,
        choices=TrainedModels.ALL,
        required=False,
        help="Sequence similarity trained model (default: %s)" % TrainedModels.MODEL_3
    )
    parser_optional.add_argument(
        "--sequence-similarity-threshold",
        dest="sequence_similarity_threshold",
        type=float,
        default=GENERATE_SEQUENCE_SIMILARITY_THRESHOLD,
        required=False,
        help="Sequence similarity threshold (default: %f). "
             "A higher threshold leads to more stringent peptide pairing." % GENERATE_SEQUENCE_SIMILARITY_THRESHOLD
    )
    parser_optional.add_argument(
        "--sequence-similarity-function",
        dest="sequence_similarity_function",
        type=str,
        default=GENERATE_SEQUENCE_SIMILARITY_FUNCTION,
        required=False,
        choices=SequenceSimilarityFunctions.ALL,
        help="Sequence similarity function (default: %s)." % GENERATE_SEQUENCE_SIMILARITY_FUNCTION
    )
    parser_optional.add_argument(
        "--random-seed",
        dest="random_seed",
        type=int,
        default=GENERATE_RANDOM_SEED,
        required=False,
        help="Random seed (default: %i)." % GENERATE_RANDOM_SEED
    )

    parser_optional_golfy = parser.add_argument_group("optional arguments (applies when '--mode golfy')")
    parser_optional_golfy.add_argument(
        "--golfy-max-iters",
        dest="golfy_max_iters",
        type=int,
        required=False,
        default=GENERATE_GOLFY_MAX_ITERS,
        help="Number of maximum iterations for golfy (default: %i)." % GENERATE_GOLFY_MAX_ITERS
    )
    parser_optional_golfy.add_argument(
        "--golfy-init-mode",
        dest="golfy_init_mode",
        type=str,
        default=GENERATE_GOLFY_INIT_MODE,
        choices=GolfyInitModes.ALL,
        required=False,
        help="Initialization mode for golfy (default: %s)." % GENERATE_GOLFY_INIT_MODE
    )

    parser_optional_sat_solver = parser.add_argument_group("optional arguments (applies when '--mode sat_solver')")
    parser_optional_sat_solver.add_argument(
        "--num-processes",
        dest="num_processes",
        type=int,
        required=False,
        default=GENERATE_NUM_PROCESSES,
        help="Number of processes (default: %i)." % GENERATE_NUM_PROCESSES
    )
    parser_optional_sat_solver.add_argument(
        "--num-peptides-per-batch",
        dest="num_peptides_per_batch",
        type=int,
        default=GENERATE_NUM_PEPTIDES_PER_BATCH,
        required=False,
        help="Number of peptides per batch (default: %i). "
             "The SAT-solver divides peptides into the specified number of peptides "
             "if the total number of peptides is bigger than the specified number "
             "(e.g. 220 peptides are divided into 2 batches of 100 peptides and 1 "
             "batch of 20 peptides if --num-peptides-per-batch is 100)." % GENERATE_NUM_PEPTIDES_PER_BATCH
    )
    parser.set_defaults(which='generate')
    return sub_parsers


def run_ace_generate_from_parsed_args(args):
    """
    Runs ACE 'generate' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   argparse.ArgumentParser object
                with the following variables:
                num_peptides
                peptides_csv_file
                num_peptides_per_pool
                num_coverage
                output_csv_file
                mode
                assign_well_ids
                plate_type
                trained_model
                sequence_similarity_function
                sequence_similarity_threshold
                random_seed
                golfy_max_iters
                golfy_init_mode
                num_processes
                num_peptides_per_batch
    """
    # Step 1. Load peptide data
    if args.peptides_csv_file is not None:
        df_peptides = pd.read_csv(args.peptides_csv_file)
        is_sequence_available = True
    else:
        data = {
            'peptide_id': [],
            'peptide_sequence': []
        }
        for i in range(1, args.num_peptides + 1):
            data['peptide_id'].append('peptide_%i' % i)
            data['peptide_sequence'].append('')
        df_peptides = pd.DataFrame(data)
        is_sequence_available = False

    # Step 2. Identify disallowed / enforced peptide pairs
    if is_sequence_available:
        # Load trained model
        trained_model_file = pkg_resources.resource_filename('acelib', 'resources/models/%s' % args.trained_model)
        ESM2_TOKENIZER = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
        ESM2_MODEL = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D", return_dict=True, output_hidden_states=True)
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        ace_eng = AceNeuralEngine(ESM2_MODEL, ESM2_TOKENIZER, device)
        ace_eng.load_weights(trained_model_file)
        preferred_peptide_pairs = ace_eng.find_paired_peptides(
            peptide_ids=df_peptides['peptide_id'].values.tolist(),
            peptide_sequences=df_peptides['peptide_sequence'].values.tolist(),
            sim_fxn=args.sequence_similarity_function,
            threshold=args.sequence_similarity_threshold
        )
        peptide_clusters = [] # todo: compute peptide clusters

        # Deduplicate
        preferred_peptide_pairs_deduped = []
        for peptide_id_1, peptide_id_2 in preferred_peptide_pairs:
            if ((peptide_id_1, peptide_id_2) not in preferred_peptide_pairs_deduped) and \
                ((peptide_id_2, peptide_id_1) not in preferred_peptide_pairs_deduped):
                preferred_peptide_pairs_deduped.append((peptide_id_1, peptide_id_2))
        logger.info('Based on our sequence similarity neural engine, here are '
                    'the peptide clusters that we will try to pool together '
                    '(%i clusters):' % len(peptide_clusters))
        for peptide_cluster in peptide_clusters:
            logger.info(peptide_cluster)
    else:
        peptide_clusters = []

    # Step 3. Generate an ELISpot configuration
    if (args.num_peptides_per_pool > 10 or args.num_coverage > 3) and \
        (args.mode == GenerateModes.SAT_SOLVER):
        logger.info('Please note that we recommend '
                    '--num-peptides-per-pool to be equal to or less than 10 and '
                    '--num-coverage to be equal to or less than 3 for --mode sat_solver. '
                    'We will still run SAT-solver but please note that this computation might '
                    'take a long time (more than 24 hours at least based on internal testing). '
                    'We recommend using --mode golfy for the requested parameters.')

    if args.mode == GenerateModes.GOLFY:
        is_valid, df_configuration = run_ace_golfy(
            df_peptides=df_peptides,
            num_peptides_per_pool=args.num_peptides_per_pool,
            num_coverage=args.num_coverage,
            random_seed=args.random_seed,
            max_iters=args.golfy_max_iters,
            init_mode=args.golfy_init_mode,
            preferred_peptide_pairs=preferred_peptide_pairs_deduped
        )
        if not is_valid:
            logger.info('The configuration generated by golfy did not result in '
                        'each peptide with a unique combination of pools. '
                        'Try running ACE again but with a higher number of '
                        '--golfy-max-iters')
    elif args.mode == GenerateModes.SAT_SOLVER:
        if len(peptide_clusters) > 0:
            df_configuration_first_coverage = ELISpot.generate_first_coverage_configuration(
                df_peptides=df_peptides,
                peptide_clusters=peptide_clusters,
                num_peptides_per_pool=args.num_peptides_per_pool
            )
            disallowed_peptide_pairs = ELISpot.compute_disallowed_peptide_pairs(
                df_configuration=df_configuration_first_coverage
            )
            is_first_coverage = False
        else:
            df_configuration_first_coverage = pd.DataFrame()
            disallowed_peptide_pairs = []
            is_first_coverage = True
        df_configuration = run_ace_sat_solver(
            df_peptides=df_peptides,
            num_peptides_per_pool=args.num_peptides_per_pool,
            num_coverage=args.num_coverage - 1,
            num_peptides_per_batch=args.num_peptides_per_batch,
            random_seed=args.random_seed,
            num_processes=args.num_processes,
            is_first_coverage=is_first_coverage,
            disallowed_peptide_pairs=disallowed_peptide_pairs
        )
        df_configuration = pd.concat([df_configuration_first_coverage, df_configuration])
    else:
        logger.error('Unknown mode: %s' % args.mode)
        exit(1)

    # Step 4. Assign plate and well IDs
    if args.assign_well_ids:
        df_configuration = ELISpot.assign_well_ids(
            df_configuration=df_configuration,
            plate_type=args.plate_type
        )

    df_configuration.to_csv(args.output_csv_file, index=False)

