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
The purpose of this python3 script is to implement utility functions.
"""


import random
import socket
from acelib.logger import get_logger


logger = get_logger(__name__)


def find_port_addr(connection):
    try:
        return connection.raddr.port
    except:
        return


def get_open_port() -> int:
    port = 1
    for i in range(10000):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            s.bind(("127.0.0.1", i))
        except socket.error as e:
             continue
        s.close()
        port = i
        break
    return port


def assign_96_well_plate_physical_ids(df_configuration, peptide_sequences):
    """
    Assigns physical plate pool IDs for a given ELIspot configuration
    for the use of 96-well plates.

    Parameters
    ----------
    df_configuration    :   DataFrame of ELIspot configuration.
    peptide_sequences   :   Peptide sequences.

    Returns
    -------
    df_configuration    :   DataFrame of ELIspot configuration with the
                            following columns added:
                            'plate_number', 'plate_well_id'
    """
    df_wells = df_configuration.loc[:,['coverage_id', 'pool_id']].drop_duplicates()
    curr_plate_number = 1
    curr_wells_used = 0
    row_prefixes = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    col_prefixes = range(1, 13)
    curr_row_idx = 0
    curr_col_idx = 0
    plate_numbers = []
    plate_well_ids = []
    for _, _ in df_wells.iterrows():
        curr_row_letter = row_prefixes[curr_row_idx]
        curr_col_num = col_prefixes[curr_col_idx]
        plate_numbers.append(curr_plate_number)
        plate_well_ids.append(str(curr_row_letter) + str(curr_col_num))
        curr_col_idx += 1
        if curr_col_idx == len(col_prefixes):
            curr_row_idx += 1
            curr_col_idx = 0
        curr_wells_used += 1
        if curr_wells_used == 96:
            curr_wells_used = 0
            curr_plate_number += 1
            curr_row_idx = 0
            curr_col_idx = 0
    df_wells['plate_number'] = plate_numbers
    df_wells['plate_well_id'] = plate_well_ids

    plate_numbers = []
    plate_well_ids = []
    peptide_sequence_values = []
    for index, row in df_configuration.iterrows():
        curr_coverage_id = row['coverage_id']
        curr_pool_id = row['pool_id']
        curr_peptide_id = row['peptide_id']
        curr_peptide_idx = int(curr_peptide_id.split('_')[1]) - 1
        df_matched = df_wells.loc[
            (df_wells['coverage_id'] == curr_coverage_id) &
            (df_wells['pool_id'] == curr_pool_id),:
        ]
        peptide_sequence_values.append(peptide_sequences[curr_peptide_idx])
        plate_numbers.append(df_matched['plate_number'].values.tolist()[0])
        plate_well_ids.append(df_matched['plate_well_id'].values.tolist()[0])

    df_configuration['sequence'] = peptide_sequence_values
    df_configuration['plate_number'] = plate_numbers
    df_configuration['plate_well_id'] = plate_well_ids
    return df_configuration
