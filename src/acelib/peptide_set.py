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
The purpose of this python3 script is to implement the PeptideSet class.
"""


import pandas as pd
from dataclasses import dataclass, field
from typing import List
from .logger import get_logger
from .peptide import Peptide


logger = get_logger(__name__)


@dataclass(frozen=True)
class PeptideSet:
    peptides: List[Peptide] = field(default_factory=list)

    def add_peptide(self, peptide: Peptide):
        """
        Add a peptide.

        Parameters:
            peptide     :   Peptide object.
        """
        self.peptides.append(peptide)

    def to_dataframe(self) -> pd.DataFrame:
        """
        Get a Pandas DataFrame.

        Returns:
            pd.DataFrame    :   A DataFrame with the following columns:

                                    - 'peptide_id'
                                    - 'peptide_sequence'
        """
        data = {
            'peptide_id': [],
            'peptide_sequence': []
        }
        for peptide in self.peptides:
            data['peptide_id'].append(peptide.id)
            data['peptide_sequence'].append(peptide.sequence)
        df = pd.DataFrame(data)
        return df
