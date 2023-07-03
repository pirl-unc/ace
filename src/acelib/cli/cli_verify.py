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
and run ACE 'verify' command.
"""


import pandas as pd
from ..logger import get_logger
from ..main import *


logger = get_logger(__name__)


def add_ace_verify_arg_parser(sub_parsers):
    """
    Adds 'verify' parser.

    Parameters
    ----------
    sub_parsers  :   argparse.ArgumentParser subparsers object.

    Returns
    -------
    argparse.ArgumentParser subparsers object.
    """
    parser = sub_parsers.add_parser(
        'verify',
        help='Verifies whether an ELISpot configuration satisfies all ACE constraints.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--configuration-csv-file",
        dest="configuration_csv_file",
        type=str,
        required=True,
        help="ELISpot configuration CSV file. "
             "Expected columns: 'coverage_id', 'pool_id', 'peptide_id'."
    )
    parser_required.add_argument(
        "--num-peptides-per-pool",
        dest="num_peptides_per_pool",
        type=int,
        required=True,
        help="Number of peptides per pool."
    )
    parser_required.add_argument(
        "--num-coverage",
        dest="num_coverage",
        type=int,
        required=True,
        help="Total coverage (i.e. number of peptide replicates)."
    )
    parser.set_defaults(which='verify')
    return sub_parsers


def run_ace_verify_from_parsed_args(args):
    """
    Runs ACE 'verify' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   argparse.ArgumentParser with the following variables:
                configuration_csv_file
                num_peptides_per_pool
                num_coverage
    """
    df_configuration = pd.read_csv(args.configuration_csv_file)
    is_optimal = run_ace_verify(
        df_configuration=df_configuration,
        num_peptides_per_pool=args.num_peptides_per_pool,
        num_coverage=args.num_coverage
    )
    if is_optimal:
        logger.info("The input ELISpot configuration meets all criteria for an "
                    "optimal configuration.")
    else:
        logger.info("The input ELISpot configuration does not meet all criteria "
                    "for an optimal configuration and is a sub-optimal configuration.")

