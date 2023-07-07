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
from ..aid_plate_reader import AIDPlateReader
from ..constants import ReadOutFileTypes
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
        help='Identify hit peptide IDs given read-outs from an ELISpot experiment.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--readout-file-type",
        dest="readout_file_type",
        type=str,
        required=True,
        help="ELISpot read-out file type (allowed options: %s)." %
             (', '.join(ReadOutFileTypes.ALL))
    )
    parser_required.add_argument(
        "--readout-files",
        dest="readout_files",
        type=str,
        action='append',
        required=True,
        help="ELISpot read-out file(s)."
             "If the readout-file-type is 'pool_id', then the expected columns of the CSV file are: 'pool_id', 'spot_count'. "
             "If the readout-file-type is 'aid_plate_reader', then the readout-files are the Excel files from the AID plate reader machine. "
             "If the readout-file type is 'aid_plate_reader' and there were pools in 2 or more plates, then supply the files in the following order: plate 1 readout file, plate 2 readout file etc."
    )
    parser_required.add_argument(
        "--configuration-csv-file",
        dest="configuration_csv_file",
        type=str,
        required=True,
        help="ELISpot configuration CSV file. "
             "Expected columns: 'pool_id', 'peptide_id', 'peptide_sequence'. Please note that if --read-out-file-type is aid_plate_reader, "
             "then the configuration must have these additional columns: 'plate_id', 'well_id'."
    )
    parser_required.add_argument(
        "--min-positive-spot-count",
        dest="min_positive_spot_count",
        type=int,
        required=True,
        help="Number of spots for a pool to be considered a positive hit."
    )
    parser_required.add_argument(
        "--output-csv-file",
        dest="output_csv_file",
        type=str,
        required=True,
        help="Output CSV file."
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
                readout_file_type
                readout_files
                configuration_csv_file
                min_positive_spot_count
                output_csv_file
    """
    # Step 1. Load the original ELISpot configuration
    df_configuration = pd.read_csv(args.configuration_csv_file)

    # Step 2. Load the read-out data
    if args.readout_file_type == ReadOutFileTypes.POOL_IDS:
        df_readout = pd.read_csv(args.readout_files[0])
        df_readout = df_readout[df_readout['spot_count'] >= args.min_positive_spot_count]
        hit_pool_ids = list(df_readout['pool_id'].unique())
    elif args.readout_file_type == ReadOutFileTypes.AID_PLATE_READER:
        if 'plate_id' not in df_configuration.columns.values.tolist():
            logger.error("The column 'plate_id' must be present in the CSV configuration file "
                         "if you supply AID plate reader XLSX read-out file(s) for deconvolution to work.")
            exit(1)
        if 'well_id' not in df_configuration.columns.values.tolist():
            logger.error("The column 'well_id' must be present in the CSV configuration file "
                         "if you supply AID plate reader XLSX read-out file(s) for deconvolution to work.")
            exit(1)

        df_hits_all = pd.DataFrame()
        plate_id = 1
        for file in args.readout_files:
            df_hits = AIDPlateReader.load_readout_file(excel_file=file, plate_id=plate_id)
            df_hits_all = pd.concat([df_hits_all, df_hits])
            plate_id += 1
        df_readout = pd.merge(df_configuration, df_hits_all, on=['plate_id', 'well_id'])
        hit_pool_ids = list(df_readout.loc[df_readout['spot_count'] >= args.min_positive_spot_count, 'pool_id'].unique())

    df_hits_max = run_ace_identify(
        hit_pool_ids=hit_pool_ids,
        df_configuration=df_configuration
    )
    df_hits_max.to_csv(args.output_csv_file, index=False)
