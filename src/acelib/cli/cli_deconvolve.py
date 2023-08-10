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


import pandas as pd
from golfy import deconvolve
from ..plate_readout import PlateReadout
from ..constants import ReadoutFileTypes
from ..deconvolution import convert_to_golfy_spotcounts, empirical_deconvolve
from ..default_parameters import *
from ..logger import get_logger
from ..main import *


logger = get_logger(__name__)


def add_ace_deconvolve_arg_parser(sub_parsers):
    """
    Adds 'deconvolve' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'deconvolve',
        help='Deconvolve hit peptide IDs given read-outs from an ELISpot experiment.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--readout-file-type",
        dest="readout_file_type",
        type=str,
        required=True,
        help="ELISpot readout file type (allowed options: %s)." %
             (', '.join(ReadoutFileTypes.ALL))
    )
    parser_required.add_argument(
        "--readout-files",
        dest="readout_files",
        type=str,
        action='append',
        required=True,
        help="ELISpot readout file(s)."
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
             "Expected columns: 'pool_id', 'peptide_id', 'peptide_sequence' "
             "in a sheet named 'block_assignment'. "
             "Expected columns: ''. "
             "Please note that if --read-out-file-type is aid_plate_reader, "
             "then the configuration must have these additional columns: 'plate_id', 'well_id'."
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
        "--mode",
        dest="mode",
        type=str,
        default=DeconvolveModes.EMPIRICAL,
        choices=DeconvolveModes.ALL,
        required=False,
        help="Deconvolution mode (default: %s)." % DeconvolveModes.EMPIRICAL
    )
    parser_optional.add_argument(
        "--min-spot-count",
        dest="min_spot_count",
        type=int,
        default=300,
        required=False,
        help="Number of spots for a pool to be considered a positive hit."
    )
    parser_optional.add_argument(
        "--min-peptide-activity",
        dest="min_peptide_activity",
        type=float,
        default=DECONVOLVE_MIN_PEPTIDE_ACTIVITY,
        required=False,
        help="Minimum estimated activity of a peptide to be considered for hit set. "
             "This value is used if when '--mode em' or '--mode lasso' (default: %f)." % DECONVOLVE_MIN_PEPTIDE_ACTIVITY
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

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser
                with the following variables:
                readout_file_type
                readout_files
                assignment_excel_file
                min_positive_spot_count
                output_excel_file
                mode
    """
    # Step 1. Load the original block assignment and design data
    block_assignment = BlockAssignment.read_excel_file(excel_file=args.assignment_excel_file)
    block_design = BlockDesign.read_excel_file(excel_file=args.assignment_excel_file)

    # Step 2. Load the readout data
    if args.readout_file_type == ReadoutFileTypes.POOL_IDS:
        df_readout = pd.read_excel(args.readout_files[0])
    elif args.readout_file_type == ReadoutFileTypes.AID_PLATE_READER:
        df_assignment = block_assignment.to_dataframe()
        if 'plate_id' not in df_assignment.columns.values.tolist():
            logger.error("The column 'plate_id' must be present in the Excel assignment file "
                         "if you supply AID plate reader XLSX readout file(s) for deconvolution.")
            exit(1)
        if 'well_id' not in df_assignment.columns.values.tolist():
            logger.error("The column 'well_id' must be present in the Excel assignment file "
                         "if you supply AID plate reader XLSX readout file(s) for deconvolution.")
            exit(1)
        plate_readouts = []
        plate_id = 1
        for file in args.readout_files:
            plate_readout = PlateReadout.read_aid_plate_reader_file(
                excel_file=file, plate_id=plate_id
            )
            plate_readouts.append(plate_readout)
            plate_id += 1
        plate_readout = PlateReadout.merge(plate_readouts=plate_readouts)
        df_readout = pd.merge(df_assignment, plate_readout.to_dataframe(), on=['plate_id', 'well_id'])

    # Step 3. Perform deconvolution
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        mode=args.mode,
        statistical_min_peptide_activity=args.min_peptide_activity,
        empirical_min_coverage=block_design.num_coverage,
        empirical_min_spot_count=args.min_spot_count,
        verbose=args.verbose
    )

    # Step 4. Write to an Excel file
    deconvolution_result.to_dataframe().to_excel(
        args.output_excel_file,
        sheet_name='deconvolved',
        index=False
    )
