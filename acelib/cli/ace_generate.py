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


from __future__ import print_function, division, absolute_import


from ..logging import get_logger
from ..default_parameters import *
from ..solver import *


logger = get_logger(__name__)


def add_ace_generate_arg_parser(sub_parsers):
    """
    Adds 'generate' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'generate',
        help='Generates an ELIspot configuration.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
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
        help="Number of peptides per pool. Please make sure this integer is a factor of the total number of peptides."
    )
    parser_required.add_argument(
        "--num_coverage",
        dest="num_coverage",
        type=int,
        required=True,
        help="Total coverage (i.e. number of peptide replicates)."
    )
    parser_required.add_argument(
        "--num_processes",
        dest="num_processes",
        type=int,
        required=True,
        default=NUM_PROCESSES,
        help="Number of processes to parallelize the computation (recommended minimum: 8)."
    )
    parser_required.add_argument(
        "--output_tsv_file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )
    parser.set_defaults(which='generate')
    return sub_parsers


def run_ace_generate_from_parsed_args(args):
    """
    Runs ACE 'generate' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser
                with the following variables:
                num_peptides
                num_peptides_per_pool
                num_coverage
                num_processes
                output_tsv_file
    """
    df_configuration = generate_assay_configuration(
        n_peptides=args.num_peptides,
        n_peptides_per_pool=args.num_peptides_per_pool,
        n_coverage=args.num_coverage,
        num_processes=args.num_processes
    )
    df_configuration.to_csv(args.output_tsv_file, sep='\t', index=False)
