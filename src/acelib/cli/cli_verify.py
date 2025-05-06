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

    Parameters:
        sub_parsers  :   argparse.ArgumentParser subparsers object.

    Returns:
        argparse.ArgumentParser subparsers object.
    """
    parser = sub_parsers.add_parser(
        'verify',
        help='Verifies whether an ELISpot assignment satisfies all ACE constraints.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--assignment-excel-file",
        dest="assignment_excel_file",
        type=str,
        required=True,
        help="ELISpot assignment Excel file. "
             "The following columns are expected to be present in a "
             "sheet named 'assignment': 'plate_id', 'well_id', 'peptide_id', 'peptide_sequence'. "
             "The following columns are expected to be present in a "
             "sheet named 'parameters': 'num_coverage', 'num_peptides_per_pool'."
    )
    parser.set_defaults(which='verify')
    return sub_parsers


def run_ace_verify_from_parsed_args(args):
    """
    Runs ACE 'verify' command using parameters from parsed arguments.

    Parameters:
        args    :   argparse.ArgumentParser with the following variables:
                    assignment_excel_file
    """
    block_assignment = BlockAssignment.read_excel_file(
        excel_file=args.assignment_excel_file
    )
    block_design = BlockDesign.read_excel_file(
        excel_file=args.assignment_excel_file
    )
    is_optimal = block_assignment.is_optimal(
        num_peptides_per_pool=block_design.num_peptides_per_pool,
        num_coverage=block_design.num_coverage
    )
    if is_optimal:
        logger.info("The input ELISpot assignment meets all criteria for an "
                    "optimal assignment.")
    else:
        logger.info("The input ELISpot assignment does not meet all criteria "
                    "for an optimal assignment and is a sub-optimal assignment.")

