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


import golfy
import pandas as pd
import random
import socket
from typing import Dict, List
from .block_assignment import BlockAssignment
from .logger import get_logger
from .peptide import Peptide


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


def convert_to_golfy_spot_counts(
        spot_counts: Dict[int,int],
        block_assignment: BlockAssignment
) -> golfy.SpotCounts:
    """
    Convert pool spot counts to golfy SpotCounts.

    Parameters:
        spot_counts             :   A mapping of pool ID to the corresponding pool spot count.
        block_assignment        :   BlockAssignment object.

    Returns:
        golfy.SpotCounts        :   golfy.SpotCounts object.
    """
    # spot_counts:
    # {
    #     coverage_id_1: {
    #         pool_id_1: spot_count_1,
    #         pool_id_2: spot_count_2,
    #         ...
    #     }
    # }
    golfy_spot_counts: Dict[int, Dict[int, int]] = {}
    for pool_id, pool in block_assignment.pools.items():
        coverage_id = pool.coverage_id - 1
        golfy_pool_id = pool_id - 1
        if coverage_id not in golfy_spot_counts:
            golfy_spot_counts[coverage_id] = {}
        if pool_id in spot_counts:
            golfy_spot_counts[coverage_id][golfy_pool_id] = spot_counts[pool_id]
    return golfy_spot_counts


def convert_dataframe_to_peptides(df_peptides: pd.DataFrame) -> List[Peptide]:
    """
    Convert pd.DataFrame to Peptides.

    Parameters:
        df_peptides :   pd.DataFrame with the following columns:
                            - 'peptide_id'
                            - 'peptide_sequence'

    Returns:
        peptides    :   List of Peptide objects.
    """
    peptides = []
    for index, value in df_peptides.iterrows():
        peptide = Peptide(id=str(value['peptide_id']), sequence=str(value['peptide_sequence']))
        peptides.append(peptide)
    return peptides


def generate_random_seed():
    return random.randint(1, 100000000)


def is_prime(num: int) -> bool:
    """
    Returns True if the supplied number is a prime number.

    Parameters:
        num         :   Number (integer).

    Returns:
        is_prime    :   True if the supplied number is a prime number.
    """
    for i in range(2, num):
        if (num % i) == 0:
            return False
    return True

