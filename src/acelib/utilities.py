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


import pandas as pd
import math
import socket
from collections import defaultdict
from typing import Dict, List, Tuple
from .logger import get_logger


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


def convert_golfy_results(golfy_assignment: Dict):
    """
    Converts golfy results into a DataFrame.

    Parameters
    ----------
    golfy_assignment    :   Dictionary.

    Returns
    -------
    df_configuration    :   DataFrame with the following columns:
                            'coverage_id'
                            'pool_id'
                            'peptide_id'
    """
    data = {
        'coverage_id': [],
        'pool_id': [],
        'peptide_id': []
    }
    curr_pool_idx = 1
    for key, value in golfy_assignment.items():
        curr_coverage_id = 'coverage_%i' % (key + 1)
        for key2, value2 in value.items():
            curr_pool_id = 'pool_%i' % curr_pool_idx
            curr_pool_idx += 1
            for p in value2:
                data['coverage_id'].append(curr_coverage_id)
                data['pool_id'].append(curr_pool_id)
                data['peptide_id'].append('peptide_%i' % (p + 1))
    return pd.DataFrame(data)


def split_peptides(
        df_peptides: pd.DataFrame,
        num_peptides_per_batch: int
) -> List[pd.DataFrame]:
    """
    Splits a DataFrame of peptides into batches.

    Parameters
    ----------
    df_peptides                 :   DataFrame with the following columns:
                                    'peptide_id'
                                    'peptide_sequence'
    num_peptides_per_batch      :   Number of peptides per batch.

    Returns
    -------
    list_df_peptides            :   List of DataFrames.
    """
    # Step 1. Calculate the number of batches
    num_batches = math.ceil(len(df_peptides) / num_peptides_per_batch)

    # Step 2. Create a list of dictionaries
    list_dict = [defaultdict(list) for i in range(0, num_batches)]

    # Step 3. Assign the peptides into batches
    for peptide_id in df_peptides['peptide_id'].unique():
        peptide_sequence = df_peptides.loc[df_peptides['peptide_id'] == peptide_id, 'peptide_sequence'].values[0]
        for i in range(0, num_batches):
            if len(list_dict[i]['peptide_id']) < num_peptides_per_batch:
                if peptide_id not in list_dict[i]['peptide_id']:
                    list_dict[i]['peptide_id'].append(peptide_id)
                    list_dict[i]['peptide_sequence'].append(peptide_sequence)
                break

    # Step 4. Convert dictionaries into DataFrames
    list_df = []
    for data_dict in list_dict:
        df_temp = pd.DataFrame(data_dict)
        list_df.append(df_temp)

    return list_df

