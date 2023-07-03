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
import random
from ortools.sat.python import cp_model
from ..constants import PlateTypes
from ..default_parameters import *
from ..logger import get_logger
from ..main import run_ace_generate
from ..utilities import convert_golfy_results, split_peptides


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
    parser_required.add_argument(
        "--num-peptides",
        dest="num_peptides",
        type=int,
        required=True,
        help="Total number of peptides."
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
        "--num-processes",
        dest="num_processes",
        type=int,
        required=True,
        default=GENERATE_NUM_PROCESSES,
        help="Number of processes (default: %i)." % GENERATE_NUM_PROCESSES
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
        "--golfy-max-iters",
        dest="golfy_max_iters",
        type=int,
        required=False,
        default=GENERATE_GOLFY_MAX_ITERS,
        help="Number of maximum iterations for golfy (default: %i)." % GENERATE_GOLFY_MAX_ITERS
    )
    parser_optional.add_argument(
        "--peptides-csv-file",
        dest="peptides_csv_file",
        type=str,
        required=False,
        help="Peptides CSV file with the following expected columns: "
             "'peptide_id', 'peptide_sequence'."
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
        "--random-seed",
        dest="random_seed",
        type=int,
        default=GENERATE_RANDOM_SEED,
        required=False,
        help="Random seed (default: %i)." % GENERATE_RANDOM_SEED
    )
    parser_optional.add_argument(
        "--num-peptides-per-batch",
        dest="num_peptides_per_batch",
        type=int,
        default=GENERATE_NUM_PEPTIDES_PER_BATCH,
        required=False,
        help="Number of peptides per batch (default: %i)." % GENERATE_NUM_PEPTIDES_PER_BATCH
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
                num_peptides_per_pool
                num_coverage
                num_processes
                output_csv_file
                peptides_csv_file
                assign_well_ids
                plate_type
                random_seed
                num_peptides_per_batch
    """
    # Step 1. Load peptide data
    if args.peptides_csv_file is not None:
        df_peptides = pd.read_csv(args.peptides_csv_file)
    else:
        data = {
            'peptide_id': [],
            'peptide_sequence': []
        }
        for i in range(1, args.num_peptides + 1):
            data['peptide_id'].append('peptide_%i' % i)
            data['peptide_sequence'].append('')
        df_peptides = pd.DataFrame(data)

    # Step 2. Identify disallowed / enforced peptide pairs
    # disallowed_peptide_pairs = dissimilarity_inference_func(df_peptides)
    # enforced_peptide_pairs = func()
    disallowed_peptide_pairs = []
    enforced_peptide_pairs = []

    # Step 3. Generate an ELISpot configuration
    status, df_configuration = run_ace_generate(
        df_peptides=df_peptides,
        num_peptides_per_pool=args.num_peptides_per_pool,
        num_coverage=args.num_coverage,
        num_processes=args.num_processes,
        random_seed=args.random_seed,
        disallowed_peptide_pairs=disallowed_peptide_pairs,
        enforced_peptide_pairs=enforced_peptide_pairs,
        assign_well_ids=args.assign_well_ids,
        plate_type=args.plate_type,
        num_peptides_per_batch=args.num_peptides_per_batch,
        golfy_max_iters=args.golfy_max_iters
    )
    df_configuration.to_csv(args.output_csv_file, index=False)

