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
and run ACE 'deconvolve' command.
"""


import os
import pandas as pd
from golfy import deconvolve
from openpyxl import load_workbook
from ..plate_readout import PlateReadout
from ..constants import ReadoutFileType
from ..deconvolution import perform_empirical_deconvolution, perform_statistical_deconvolution
from ..defaults import *
from ..logger import get_logger
from ..main import *
from ..utilities import convert_to_golfy_spot_counts


logger = get_logger(__name__)


def add_ace_deconvolve_arg_parser(sub_parsers):
    """
    Adds 'deconvolve' parser.

    Parameters:
        sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns:
        An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'deconvolve',
        help='Deconvolve hit peptide IDs given read-outs from a pooled ELISpot experiment.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--readout-file-type",
        dest="readout_file_type",
        type=str,
        required=True,
        choices=[str(f) for f in ReadoutFileType],
        help="ELISpot readout file type (allowed options: %s)." %
             (', '.join([str(f) for f in ReadoutFileType]))
    )
    parser_required.add_argument(
        "--readout-files",
        dest="readout_files",
        type=str,
        action='append',
        required=True,
        help="ELISpot readout file(s). "
             "If the readout-file-type is 'pool_id', then the expected columns of the Excel file are: 'pool_id', 'spot_count'. "
             "If the readout-file-type is 'aid_plate_reader', then the readout-files are the Excel files from the AID plate reader machine. "
             "If the readout-file type is 'aid_plate_reader' and there were pools in 2 or more plates, then supply the files in the following order: plate 1 readout file, plate 2 readout file etc."
    )
    parser_required.add_argument(
        "--assignment-excel-file",
        dest="assignment_excel_file",
        type=str,
        required=True,
        help="ELISpot assignment Excel file. "
             "Expected columns: 'peptide_id', 'peptide_sequence', 'plate_id', 'well_id' "
             "in a sheet named 'block_assignment'. "
    )
    parser_required.add_argument(
        "--min-pool-spot-count",
        dest="min_pool_spot_count",
        type=int,
        required=True,
        help="Minimum spot count for a pool to be considered a positive pool."
    )
    parser_required.add_argument(
        "--output-excel-file",
        dest="output_excel_file",
        type=str,
        required=True,
        help="Output deconvolution Excel file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--method",
        dest="method",
        type=str,
        default=DEFAULT_DECONVOLVE_METHOD,
        choices=[str(m) for m in DeconvolutionMethod],
        required=False,
        help="Deconvolution method (default: %s)." % DEFAULT_DECONVOLVE_METHOD
    )
    parser_optional.add_argument(
        "--verbose",
        dest="verbose",
        type=bool,
        required=False,
        default=True,
        help="If True, prints messages. Otherwise, messages are not printed (default: True)."
    )
    parser.set_defaults(which='deconvolve')
    return sub_parsers


def run_ace_deconvolve_from_parsed_args(args):
    """
    Runs ACE 'deconvolve' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:

                        - readout_file_type
                        - readout_files
                        - assignment_excel_file
                        - min_positive_pool_spot_count
                        - output_excel_file
                        - statistical_deconvolution_method
    """
    # Step 1. Load the original block assignment and design data
    block_assignment = BlockAssignment.read_excel_file(excel_file=args.assignment_excel_file)
    block_design = BlockDesign.read_excel_file(excel_file=args.assignment_excel_file)

    # Step 2. Load the readout data
    if ReadoutFileType(args.readout_file_type) == ReadoutFileType.POOL_IDS:
        plate_readout = PlateReadout.read_pool_id_file(file=args.readout_files[0])
    elif ReadoutFileType(args.readout_file_type) == ReadoutFileType.AID_PLATE_READER:
        plate_readouts = []
        plate_id = 1
        for file in args.readout_files:
            plate_readout = PlateReadout.read_aid_plate_reader_file(
                excel_file=file, plate_id=plate_id
            )
            plate_readouts.append(plate_readout)
            plate_id += 1
        plate_readout = PlateReadout.merge(plate_readouts=plate_readouts)
    else:
        raise Exception('Unsupported read-out file type: %s' % args.readout_file_type)

    # Step 3. Get pool IDs
    plate_readout.assign_pool_ids(block_assignment=block_assignment)
    df_readout = plate_readout.to_dataframe()

    # Step 4. Perform deconvolution
    deconvolved_peptide_set = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        method=DeconvolutionMethod(args.method),
        min_coverage=block_design.num_coverage,
        min_pool_spot_count=args.min_pool_spot_count,
        verbose=args.verbose
    )

    # Step 5. Write to an Excel file
    df_deconvolution = deconvolved_peptide_set.to_dataframe()
    df_metadata = deconvolved_peptide_set.metadata_dataframe()

    with pd.ExcelWriter(args.output_excel_file, engine='openpyxl', mode='w') as writer:
        df_deconvolution.to_excel(
            writer,
            sheet_name='deconvolution_results',
            index=False
        )
        df_metadata.to_excel(
            writer,
            sheet_name='metadata',
            index=False
        )

