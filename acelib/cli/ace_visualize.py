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
and run ACE 'visualize' command.
"""


from __future__ import print_function, division, absolute_import


import pandas as pd
from ..logging import get_logger
from ..visualization import *


logger = get_logger(__name__)


def add_ace_visualize_arg_parser(sub_parsers):
    """
    Adds 'visualize' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'visualize',
        help='Visualizes an ELIspot configuration.'
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
        "--output_pdf_file",
        dest="output_pdf_file",
        type=str,
        required=True,
        help="Output PDF file."
    )
    parser.set_defaults(which='visualize')
    return sub_parsers


def run_ace_visualize_from_parsed_args(args):
    """
    Runs ACE 'visualize' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser
                with the following variables:
                configuration_tsv_file
                output_pdf_file
    """
    df_configuration = pd.read_csv(args.configuration_tsv_file, sep='\t')
    plt = plot_configuration_table(
        df_configuration=df_configuration,
        save_figure=True,
        output_pdf_file=args.output_pdf_file
    )
    plt.savefig(args.output_pdf_file)

