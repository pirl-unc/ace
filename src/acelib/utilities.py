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
import random
import socket
from collections import defaultdict
from typing import Dict, List, Tuple
from .logger import get_logger
from .types import *


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


def convert_peptides_to_dataframe(peptides: Peptides) -> pd.DataFrame:
    """
    Convert Peptides to pd.DataFrame.

    Parameters
    ----------
    peptides    :   Peptides.

    Returns
    -------
    df_peptides :   pd.DataFrame with the following columns:
                    'peptide_id'
                    'peptide_sequence'
    """
    data = {
        'peptide_id': [],
        'peptide_sequence': []
    }
    for peptide_id, peptide_sequence in peptides:
        data['peptide_id'].append(peptide_id)
        data['peptide_sequence'].append(peptide_sequence)
    return pd.DataFrame(data)


def convert_dataframe_to_peptides(df_peptides: pd.DataFrame) -> Peptides:
    """
    Convert pd.DataFrame to Peptides.

    Parameters
    ----------
    df_peptides :   pd.DataFrame with the following columns:
                    'peptide_id'
                    'peptide_sequence'

    Returns
    -------
    peptides    :   Peptides.
    """
    peptides = []
    for index, value in df_peptides.iterrows():
        peptides.append((value['peptide_id'], value['peptide_sequence']))
    return peptides


def generate_random_seed():
    return random.randint(1, 100000000)


def is_prime(num: int) -> bool:
    """
    Returns True if the supplied number is a prime number.

    Parameters
    ----------
    num         :   Number (integer).

    Returns
    -------
    is_prime    :   True if the supplied number is a prime number.
    """
    for i in range(2, num):
        if (num % i) == 0:
            return False
    return True

