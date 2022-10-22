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


from __future__ import print_function, division, absolute_import


import pandas as pd
from ..logger import get_logger
from ..main import *


logger = get_logger(__name__)


def add_ace_verify_arg_parser(sub_parsers):
    """
    Adds 'verify' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'verify',
        help='Verifies whether an ELIspot configuration satisfies all ACE constraints.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--configuration_tsv_file",
        dest="configuration_tsv_file",
        type=str,
        required=True,
        help="ELIspot configuration TSV file. "
             "Expected columns: 'coverage_id', 'pool_id', 'peptide_id'."
    )
    parser_required.add_argument(
        "--num_peptides",
        dest="num_peptides",
        type=int,
        required=True,
        help="Total number of peptides."
    )
    parser_required.add_argument(
        "--num_peptides_per_pool",
        dest="num_peptides_per_pool",
        type=int,
        required=True,
        help="Number of peptides per pool."
    )
    parser_required.add_argument(
        "--num_coverage",
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
    args    :   An instance of argparse.ArgumentParser
                with the following variables:
                configuration_tsv_file
                num_peptides
                num_peptides_per_pool
                num_coverage
    """
    df_configuration = pd.read_csv(args.configuration_tsv_file, sep='\t')
    is_valid = run_ace_verify(
        df_configuration=df_configuration,
        n_peptides=args.num_peptides,
        n_peptides_per_pool=args.num_peptides_per_pool,
        n_coverage=args.num_coverage
    )
    if is_valid:
        logger.info("ELIspot configuration meets all ACE constraints and is valid.")
    else:
        logger.info("ELIspot configuration does not meet all ACE constraints and is not valid.")

