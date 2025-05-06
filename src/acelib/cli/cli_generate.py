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
from random import choices

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
from ..defaults import *
from ..logger import get_logger
from ..main import run_ace_sat_solver, run_ace_golfy, run_ace_generate
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
        "--peptides-file",
        dest="peptides_file",
        type=str,
        help="Peptides CSV or Excel file with the following columns: "
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
        help="Output Excel (.xlsx) file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--num-plate-wells",
        dest="num_plate_wells",
        type=int,
        default=DEFAULT_GENERATE_NUM_PLATE_WELLS,
        choices=[int(n.value) for n in NumPlateWells],
        required=False,
        help="Number of wells on plate. Allowed values: %s (default: %i)." %
             (', '.join([str(n.value) for n in NumPlateWells]), int(NumPlateWells.WELLS_96.value))
    )
    parser_optional.add_argument(
        "--mode",
        dest="mode",
        type=str,
        default=DEFAULT_GENERATE_MODE,
        choices=[str(m) for m in GenerateMode],
        required=False,
        help="Configuration generation mode. Allowed values: %s (default: %s)." %
             (', '.join([str(m) for m in GenerateMode]), GenerateMode.GOLFY)
    )
    parser_optional.add_argument(
        "--cluster-peptides",
        dest="cluster_peptides",
        type=bool,
        default=DEFAULT_GENERATE_CLUSTER_PEPTIDES,
        required=False,
        help="Cluster peptides if set to true (default: %r)." % DEFAULT_GENERATE_CLUSTER_PEPTIDES
    )
    parser_optional.add_argument(
        "--sequence-similarity-function",
        dest="sequence_similarity_function",
        type=str,
        default=DEFAULT_GENERATE_SEQUENCE_SIMILARITY_FUNCTION,
        required=False,
        choices=[str(f) for f in SequenceSimilarityFunction],
        help="Sequence similarity function. Allowed values: %s (default: %s)." %
             (', '.join([str(f) for f in SequenceSimilarityFunction]), DEFAULT_GENERATE_SEQUENCE_SIMILARITY_FUNCTION)
    )
    parser_optional.add_argument(
        "--sequence-similarity-threshold",
        dest="sequence_similarity_threshold",
        type=float,
        default=DEFAULT_GENERATE_SEQUENCE_SIMILARITY_THRESHOLD,
        required=False,
        help="Sequence similarity threshold (default: %f). "
             "A higher threshold leads to more stringent peptide pairing." % DEFAULT_GENERATE_SEQUENCE_SIMILARITY_THRESHOLD
    )
    # Golfy optional parameters
    parser_optional_golfy = parser.add_argument_group("optional arguments (applies when '--mode golfy')")
    parser_optional_golfy.add_argument(
        "--golfy-random-seed",
        dest="golfy_random_seed",
        type=int,
        default=DEFAULT_GENERATE_GOLFY_RANDOM_SEED,
        required=False,
        help="Random seed for golfy (default: %i)." % DEFAULT_GENERATE_GOLFY_RANDOM_SEED
    )
    parser_optional_golfy.add_argument(
        "--golfy-max-iters",
        dest="golfy_max_iters",
        type=int,
        required=False,
        default=DEFAULT_GENERATE_GOLFY_MAX_ITERS,
        help="Number of maximum iterations for golfy (default: %i)." % DEFAULT_GENERATE_GOLFY_MAX_ITERS
    )
    parser_optional_golfy.add_argument(
        "--golfy-strategy",
        dest="golfy_strategy",
        type=str,
        default=DEFAULT_GENERATE_GOLFY_STRATEGY,
        choices=[str(s) for s in GolfyStrategy],
        required=False,
        help="Strategy for golfy. Allowed values: %s (default: %s)." %
             ([str(s) for s in GolfyStrategy], DEFAULT_GENERATE_GOLFY_STRATEGY)
    )
    parser_optional_golfy.add_argument(
        "--golfy-allow-extra-pools",
        dest="golfy_allow_extra_pools",
        type=eval,
        default=DEFAULT_GENERATE_GOLFY_ALLOW_EXTRA_POOLS,
        choices=[True, False],
        required=False,
        help="Allow extra pools for golfy (default: %r)." % DEFAULT_GENERATE_GOLFY_ALLOW_EXTRA_POOLS
    )

    # CP-SAT solver optional parameters
    parser_optional_sat_solver = parser.add_argument_group("optional arguments (applies when '--mode cpsat_solver')")
    parser_optional_sat_solver.add_argument(
        "--cpsat-solver-num-processes",
        dest="cpsat_solver_num_processes",
        type=int,
        required=False,
        default=DEFAULT_GENERATE_CPSAT_SOLVER_NUM_PROCESSES,
        help="Number of processes for CP-SAT solver (default: %i)." % DEFAULT_GENERATE_CPSAT_SOLVER_NUM_PROCESSES
    )
    parser_optional_sat_solver.add_argument(
        "--cpsat-solver-shuffle-iters",
        dest="cpsat_solver_shuffle_iters",
        type=int,
        required=False,
        default=DEFAULT_GENERATE_CPSAT_SOLVER_SHUFFLE_ITERS,
        help="Number of iterations to shuffle pool IDs to minimize "
             "number of non-unique pool assignment violations for CP-SAT solver (default: %i)." % DEFAULT_GENERATE_CPSAT_SOLVER_SHUFFLE_ITERS
    )
    parser_optional_sat_solver.add_argument(
        "--cpsat-solver-max-peptides-per-block",
        dest="cpsat_solver_max_peptides_per_block",
        type=int,
        default=DEFAULT_GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_BLOCK,
        required=False,
        help="Maximum number of peptides per block (default: %i). "
             "The CPSAT-solver divides peptides into the specified number of peptides "
             "if the total number of peptides is bigger than the specified number "
             "(e.g. 220 peptides are divided into 2 blocks of 100 peptides and 1 "
             "block of 20 peptides if --max-peptides-per-block is 100). "
             "Increasing this number from the current default value will likely "
             "make the computation intractable so it is recommended that you "
             "keep this at %i." %
             (DEFAULT_GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_BLOCK, DEFAULT_GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_BLOCK)
    )
    parser_optional_sat_solver.add_argument(
        "--cpsat-solver-max-peptides-per-pool",
        dest="cpsat_solver_max_peptides_per_pool",
        type=int,
        default=DEFAULT_GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_POOL,
        required=False,
        help="Maximum number of peptides per pool (default: %i). "
             "Increasing this number from the current default value will likely "
             "make the computation intractable so it is recommended that you "
             "keep this at %i." %
             (DEFAULT_GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_POOL, DEFAULT_GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_POOL)
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
                peptides_file
                design_csv_file
                num_peptides_per_pool
                num_coverage
                output_excel_file
                mode
                plate_size
                sequence_similarity_function
                sequence_similarity_threshold
                golfy_random_seed
                golfy_max_iters
                golfy_strategy
                golfy_allow_extra_pools
                cpsat_solver_num_processes
                cpsat_solver_shuffle_iters
                cpsat_solver_max_peptides_per_block
                cpsat_solver_max_peptides_per_pool
                verbose
    """
    # Step 1. Load peptide data
    if args.peptides_file is not None:
        file_name, file_extension = os.path.splitext(args.peptides_file)
        if file_extension == '.xlsx':
            peptides = convert_dataframe_to_peptides(df_peptides=pd.read_excel(args.peptides_file))
        elif file_extension == '.csv':
            peptides = convert_dataframe_to_peptides(df_peptides=pd.read_csv(args.peptides_file))
        else:
            logger.error("--peptides-file must be either a .CSV or .XLSX file. Expected headers: 'peptide_id' and 'peptide_sequence'.")
            exit(1)
        cluster_peptides = args.cluster_peptides
    else:
        peptides = []
        for i in range(1, args.num_peptides + 1):
            peptides.append(('peptide_%i' % i, ''))
        cluster_peptides = False

    # Step 2. Check input parameters
    if len(peptides) % args.num_peptides_per_pool != 0 and args.mode == GenerateMode.CPSAT_SOLVER:
        logger.error('The number of peptides per pool must be divisible by the total number of peptides to for CP-SAT solver.')
        exit(1)  

    # Step 3. Generate block design
    block_assignment, block_design = run_ace_generate(
        peptides=peptides,
        num_peptides_per_pool=args.num_peptides_per_pool,
        num_coverage=args.num_coverage,
        trained_model_file=str(resources.path('acelib.resources.models', 'trained_model_w_data_augmentation_b3000.pt')),
        cluster_peptides=cluster_peptides,
        mode=GenerateMode(args.mode),
        sequence_similarity_function=SequenceSimilarityFunction(args.sequence_similarity_function),
        sequence_similarity_threshold=args.sequence_similarity_threshold,
        golfy_random_seed=args.golfy_random_seed,
        golfy_strategy=GolfyStrategy(args.golfy_strategy),
        golfy_max_iters=args.golfy_max_iters,
        golfy_allow_extra_pools=args.golfy_allow_extra_pools,
        cpsat_solver_num_processes=args.cpsat_solver_num_processes,
        cpsat_solver_shuffle_iters=args.cpsat_solver_shuffle_iters,
        cpsat_solver_max_peptides_per_block=args.cpsat_solver_max_peptides_per_block,
        cpsat_solver_max_peptides_per_pool=args.cpsat_solver_max_peptides_per_pool,
        num_plate_wells=NumPlateWells(args.num_plate_wells),
        verbose=args.verbose
    )

    # Step 4. Write design and assignment to an Excel file
    df_assignment = block_assignment.to_dataframe()
    df_assignment.drop('pool_id', axis=1, inplace=True)
    df_assignment = df_assignment[['peptide_id', 'peptide_sequence', 'plate_id', 'well_id', 'coverage_id']]
    df_assignment_bench_ready = block_assignment.to_bench_ready_dataframe()
    with pd.ExcelWriter(args.output_excel_file, engine='openpyxl') as writer:
        df_assignment.to_excel(writer, sheet_name='assignment', index=False)
        df_assignment_bench_ready.to_excel(writer, sheet_name='assignment_bench_ready', index=False)
        block_design.peptides_dataframe.to_excel(writer, sheet_name='peptides', index=False)
        block_design.preferred_peptide_pairs_dataframe.to_excel(writer, sheet_name='preferred_peptide_pairs', index=False)
        block_design.metadata_dataframe.to_excel(writer, sheet_name='parameters', index=False)

