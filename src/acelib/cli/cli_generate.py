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


import argparse
import math
import pandas as pd
import os
import random
import torch
from importlib import resources
from pathlib import Path
from ortools.sat.python import cp_model
from transformers import BertModel, BertTokenizer
from transformers import AutoTokenizer, AutoModelForMaskedLM
from ..block_assignment import BlockAssignment
from ..block_design import BlockDesign
from ..constants import *
from ..default_parameters import *
from ..logger import get_logger
from ..main import run_ace_sat_solver, run_ace_golfy
from ..utilities import *
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
        "--peptides-excel-file",
        dest="peptides_excel_file",
        type=str,
        help="Peptides Excel file with the following columns: "
             "'peptide_id', 'peptide_sequence'. Please note that only "
             "either this parameter or '--num-peptides' can be supplied."
    )
    parser_required.add_argument(
        "--num-peptides-per-pool",
        dest="num_peptides_per_pool",
        type=int,
        required=True,
        help="Number of peptides per pool (i.e. well). "
             "Please make sure this number is divisible by the total number of peptides."
    )
    parser_required.add_argument(
        "--num-coverage",
        dest="num_coverage",
        type=int,
        required=True,
        help="Total coverage (i.e. number of peptide replicates)."
    )
    parser_required.add_argument(
        "--output-excel-file",
        dest="output_excel_file",
        type=str,
        required=True,
        help="Output assignment Excel file."
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
    parser_optional_golfy.add_argument(
        "--golfy-allow-extra-pools",
        dest="golfy_allow_extra_pools",
        type=eval,
        default=GENERATE_GOLFY_ALLOW_EXTRA_POOLS,
        choices=[True, False],
        required=False,
        help="Allow extra pools for golfy (default: %r)." % GENERATE_GOLFY_ALLOW_EXTRA_POOLS
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
        "--shuffle-iters",
        dest="shuffle_iters",
        type=int,
        required=False,
        default=GENERATE_SHUFFLE_ITERS,
        help="Number of iterations to shuffle pool IDs to minimize "
             "number of non-unique pool assignment violations (default: %i)." % GENERATE_SHUFFLE_ITERS
    )
    parser_optional_sat_solver.add_argument(
        "--max-peptides-per-block",
        dest="max_peptides_per_block",
        type=int,
        default=GENERATE_MAX_PEPTIDES_PER_BLOCK,
        required=False,
        help="Maximum number of peptides per block (default: %i). "
             "This parameter applies when --mode sat_solver. "
             "The SAT-solver divides peptides into the specified number of peptides "
             "if the total number of peptides is bigger than the specified number "
             "(e.g. 220 peptides are divided into 2 blocks of 100 peptides and 1 "
             "block of 20 peptides if --max-peptides-per-block is 100). "
             "Increasing this number from the current default value will likely "
             "make the computation intractable so it is recommended that you "
             "keep this at %i." %
             (GENERATE_MAX_PEPTIDES_PER_BLOCK, GENERATE_MAX_PEPTIDES_PER_BLOCK)
    )
    parser_optional_sat_solver.add_argument(
        "--max-peptides-per-pool",
        dest="max_peptides_per_pool",
        type=int,
        default=GENERATE_MAX_PEPTIDES_PER_POOL,
        required=False,
        help="Maximum number of peptides per pool (default: %i). "
             "This parameter applies when --mode sat_solver. "
             "Increasing this number from the current default value will likely "
             "make the computation intractable so it is recommended that you "
             "keep this at %i." %
             (GENERATE_MAX_PEPTIDES_PER_POOL, GENERATE_MAX_PEPTIDES_PER_POOL)
    )
    parser_optional_sat_solver.add_argument(
        "--verbose",
        dest="verbose",
        type=bool,
        required=False,
        default=True,
        help="If True, prints messages. Otherwise, messages are not printed (default: True)."
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
                peptides_excel_file
                design_csv_file
                num_peptides_per_pool
                num_coverage
                output_excel_file
                mode
                assign_well_ids
                plate_type
                sequence_similarity_function
                sequence_similarity_threshold
                random_seed
                golfy_max_iters
                golfy_init_mode
                golfy_allow_extra_pools
                num_processes
                shuffle_iters
                max_peptides_per_block
                max_peptides_per_pool
                verbose
    """
    # Step 1. Load peptide data
    if args.peptides_excel_file is not None:
        peptides = convert_dataframe_to_peptides(df_peptides=pd.read_excel(args.peptides_excel_file))
        is_sequence_available = True
    else:
        peptides = []
        for i in range(1, args.num_peptides + 1):
            peptides.append(('peptide_%i' % i, ''))
        is_sequence_available = False

    # Step 2. Check input parameters
    if len(peptides) % args.num_peptides_per_pool != 0:
        logger.error('The number of peptides per pool must be divisible by the total number of peptides.')
        exit(1)

    # Step 3. Identify pairs of similar peptides
    if is_sequence_available:
        # Load trained model
        trained_model_file = resources.path('acelib.resources.models', 'seq_sim_trained_model.pt')
        ESM2_TOKENIZER = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
        ESM2_MODEL = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D", return_dict=True, output_hidden_states=True)
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        ace_eng = AceNeuralEngine(ESM2_MODEL, ESM2_TOKENIZER, device)
        ace_eng.load_weights(trained_model_file)

        # Identify pairs of similar peptides
        preferred_peptide_pairs = ace_eng.find_paired_peptides(
            peptide_ids=[p[0] for p in peptides],
            peptide_sequences=[p[1] for p in peptides],
            sim_fxn=args.sequence_similarity_function,
            threshold=args.sequence_similarity_threshold
        )
        if args.verbose:
            logger.info('Based on our sequence similarity neural engine, here are '
                        'the peptide pairs that we will try to pool together '
                        '(%i pairs):' % len(preferred_peptide_pairs))
            logger.info('peptide ID, peptide ID: similarity score')
            for peptide_id_1, peptide_id_2, score in preferred_peptide_pairs:
                logger.info('%s, %s: %f' % (peptide_id_1, peptide_id_2, score))
    else:
        preferred_peptide_pairs = []
    preferred_peptide_pairs = [(p1, p2) for p1, p2, score in preferred_peptide_pairs]

    # Step 4. Generate a block design
    block_design = BlockDesign(
        peptides=peptides,
        num_peptides_per_pool=args.num_peptides_per_pool,
        num_coverage=args.num_coverage,
        max_peptides_per_block=args.max_peptides_per_block,
        disallowed_peptide_pairs=[],
        preferred_peptide_pairs=preferred_peptide_pairs
    )

    # Step 5. Generate a block assignment
    if args.mode == GenerateModes.GOLFY:
        block_assignment = run_ace_golfy(
            block_design=block_design,
            random_seed=args.random_seed,
            max_iters=args.golfy_max_iters,
            init_mode=args.golfy_init_mode,
            allow_extra_pools=args.golfy_allow_extra_pools,
            verbose=args.verbose
        )
    elif args.mode == GenerateModes.SAT_SOLVER:
        block_assignment = run_ace_sat_solver(
            block_design=block_design,
            max_peptides_per_pool=args.max_peptides_per_pool,
            num_processes=args.num_processes,
            shuffle_iters=args.shuffle_iters,
            verbose=args.verbose
        )
    else:
        logger.error('Unknown mode: %s' % args.mode)
        exit(1)

    # Step 6. Check if the block assignment is optimal
    # re-assign value because it could have been decremented
    # to account for preferred peptide neighbors (1x coverage).
    block_design.num_coverage = args.num_coverage
    block_assignment.is_optimal(
        num_coverage=block_design.num_coverage,
        num_peptides_per_pool=block_design.num_peptides_per_pool,
        verbose=args.verbose
    )

    # Step 7. Assign plate and well IDs
    if args.assign_well_ids:
        df_assignment = BlockAssignment.assign_well_ids(
            df_assignment=block_assignment.to_dataframe(),
            plate_type=args.plate_type
        )
    else:
        df_assignment = block_assignment.to_dataframe()

    # Step 8. Write design and assignment to an Excel file
    df_design = block_design.to_dataframe()
    writer = pd.ExcelWriter(args.output_excel_file, engine='openpyxl')
    df_assignment.to_excel(writer, sheet_name='block_assignment', index=False)
    df_design.to_excel(writer, sheet_name='block_design', index=False)
    writer.save()
