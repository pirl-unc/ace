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
        default=DeconvolutionMethods.CONSTRAINED_EM,
        choices=DeconvolutionMethods.ALL,
        required=False,
        help="Deconvolution method (default: %s)." % DeconvolutionMethods.CONSTRAINED_EM
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
                min_positive_pool_spot_count
                output_excel_file
                statistical_deconvolution_method
    """
    # Step 1. Load the original block assignment and design data
    block_assignment = BlockAssignment.read_excel_file(excel_file=args.assignment_excel_file)
    block_design = BlockDesign.read_excel_file(excel_file=args.assignment_excel_file)

    # Step 2. Load the readout data
    df_assignment = block_assignment.to_dataframe()
    if args.readout_file_type == ReadoutFileTypes.POOL_IDS:
        file_name, file_extension = os.path.splitext(args.readout_files[0])
        if file_extension == '.xlsx':
            df_readout = pd.read_excel(args.readout_files[0])
        elif file_extension == '.csv':
            df_readout = pd.read_csv(args.readout_files[0])
        else:
            logger.error("--readout-files must be either a .CSV or .XLSX file. Expected headers: 'plate_id', 'well_id', and 'spot_count'.")
            exit(1)
    elif args.readout_file_type == ReadoutFileTypes.AID_PLATE_READER:
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

    # Step 3. Get pool IDs
    pool_ids = []
    well_ids_dict = {}  # key   = <pool_id>
                        # value = <plate_id>-<well_id>
    for index, row in df_readout.iterrows():
        curr_plate_id = row['plate_id']
        curr_well_id = row['well_id']
        curr_pool_id = df_assignment.loc[
            (df_assignment['plate_id'] == curr_plate_id) &
            (df_assignment['well_id'] == curr_well_id),'pool_id'].values[0]
        pool_ids.append(curr_pool_id)
        well_ids_dict[curr_pool_id] = '%s-%s' % (curr_plate_id, curr_well_id)
    df_readout['pool_id'] = pool_ids

    # Step 4. Perform deconvolution
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        method=args.method,
        min_coverage=block_design.num_coverage,
        min_pool_spot_count=args.min_pool_spot_count,
        verbose=args.verbose
    )

    # Step 5. Convert pool IDs to well IDs
    data = {
        'peptide_id': [],
        'peptide_sequence': [],
        'estimated_peptide_spot_count': [],
        'hit_well_ids': [],
        'hit_well_ids_count': [],
        'deconvolution_result': []
    }
    for index, row in deconvolution_result.to_dataframe().iterrows():
        well_ids = []
        for curr_pool_id in row['hit_pool_ids'].split(';'):
            if curr_pool_id != '':
                well_ids.append(well_ids_dict[int(curr_pool_id)])
        data['peptide_id'].append(row['peptide_id'])
        data['peptide_sequence'].append(row['peptide_sequence'])
        data['estimated_peptide_spot_count'].append(row['estimated_peptide_spot_count'])
        data['hit_well_ids'].append(','.join(well_ids))
        data['hit_well_ids_count'].append(len(well_ids))
        data['deconvolution_result'].append(row['deconvolution_result'])
    df_results = pd.DataFrame(data)

    # Step 5. Write to an Excel file
    df_results.to_excel(
        args.output_excel_file,
        sheet_name='deconvolution_results',
        index=False
    )
