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
and run ACE 'identify' command.
"""


from __future__ import print_function, division, absolute_import


import pandas as pd
from ..default_parameters import *
from ..logger import get_logger
from ..main import *


logger = get_logger(__name__)


def add_ace_identify_arg_parser(sub_parsers):
    """
    Adds 'identify' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'identify',
        help='Identify hit peptide IDs given read-outs from an ELIspot experiment.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--readout_tsv_file",
        dest="readout_tsv_file",
        type=str,
        required=True,
        help="ELIspot read-out TSV file. "
             "Expected columns: 'pool_id', 'spot_count'."
    )
    parser_required.add_argument(
        "--configuration_tsv_file",
        dest="configuration_tsv_file",
        type=str,
        required=True,
        help="ELIspot configuration TSV file. "
             "Expected columns: 'coverage_id', 'pool_id', 'peptide_id'."
    )
    parser_required.add_argument(
        "--min_positive_spot_count",
        dest="min_positive_spot_count",
        type=int,
        required=True,
        default=MIN_POSITIVE_SPOT_COUNT,
        help="Number of spots for a pool to be considered a hit."
    )
    parser_required.add_argument(
        "--output_tsv_file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )
    parser.set_defaults(which='identify')
    return sub_parsers


def run_ace_identify_from_parsed_args(args):
    """
    Runs ACE 'identify' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser
                with the following variables:
                readout_tsv_file
                configuration_tsv_file
                min_positive_spot_count
                output_tsv_file
    """
    df_readout = pd.read_csv(args.readout_tsv_file)
    df_configuration = pd.read_csv(args.configuration_tsv_file)
    df_hits = run_ace_identify(
        df_readout=df_readout,
        df_configuration=df_configuration,
        min_positive_spot_count=args.min_positive_spot_count
    )
    df_hits.to_csv(args.output_tsv_file, sep='\t', index=False)
