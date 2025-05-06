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
The purpose of this python3 script is to implement the HitPeptideSet class.
"""


import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, Mapping, List, Optional
from .constants import DeconvolutionResult, DeconvolutionMethod
from .deconvolved_peptide import DeconvolvedPeptide
from .logger import get_logger
from .plate_well import PlateWell


logger = get_logger(__name__)


@dataclass
class DeconvolvedPeptideSet:
    min_pool_spot_count: Optional[float] = None
    min_coverage: Optional[int] = None
    min_peptide_spot_count: Optional[float] = None
    background_spot_count: Optional[float] = None
    pool_spot_counts: Optional[Dict[int,float]] = None # key = pool ID, value = pool spot count
    deconvolution_method: Optional[DeconvolutionMethod] = None
    plate_map: Optional[Mapping[int, PlateWell]] = field(
        default_factory=dict,
        metadata={"doc": "Mapping from a pool ID to plate and well IDs."}
    )
    deconvolved_peptides: List[DeconvolvedPeptide] = field(default_factory=list)

    def add(self, deconvolved_peptide: DeconvolvedPeptide):
        """
        Add a DeconvolvedPeptide object.

        Parameters:
            deconvolved_peptide     :   DeconvolvedPeptide object.
        """
        self.deconvolved_peptides.append(deconvolved_peptide)

    def get_confident_hits(self) -> List[DeconvolvedPeptide]:
        """Return peptides with result == CONFIDENT_HIT."""
        return [p for p in self.deconvolved_peptides if p.result == DeconvolutionResult.CONFIDENT_HIT]

    def get_candidate_hits(self) -> List[DeconvolvedPeptide]:
        """Return peptides with result == CANDIDATE_HIT."""
        return [p for p in self.deconvolved_peptides if p.result == DeconvolutionResult.CANDIDATE_HIT]

    def get_non_hits(self) -> List[DeconvolvedPeptide]:
        """Return peptides with result == NOT_A_HIT."""
        return [p for p in self.deconvolved_peptides if p.result == DeconvolutionResult.NOT_A_HIT]

    def metadata_dataframe(self) -> pd.DataFrame:
        data = {
            'min_pool_spot_count': [self.min_pool_spot_count],
            'min_coverage': [self.min_coverage],
            'min_peptide_spot_count': [self.min_peptide_spot_count],
            'background_spot_count': [self.background_spot_count],
            'deconvolution_method': [self.deconvolution_method]
        }
        return pd.DataFrame(data)

    def to_dataframe(self) -> pd.DataFrame:
        """
        Get a Pandas DataFrame.

        Returns:
            pd.DataFrame    : A DataFrame with the following columns:

                                - 'peptide_id'
                                - 'peptide_sequence'
                                - 'estimated_peptide_spot_count'
                                - 'hit_pool_ids'
                                - 'hit_pools_count'
                                - 'hit_plate_well_ids',
                                - 'hit_pool_spot_counts',
                                - 'deconvolution_result'
        """
        data = {
            'peptide_id': [],
            'peptide_sequence': [],
            'estimated_peptide_spot_count': [],
            'hit_pool_ids': [],
            'hit_pools_count': [],
            'hit_plate_well_ids': [],
            'hit_pool_spot_counts': [],
            'deconvolution_result': []
        }

        for deconvolved_peptide in self.deconvolved_peptides:
            data['peptide_id'].append(deconvolved_peptide.id)
            data['peptide_sequence'].append(deconvolved_peptide.sequence)
            data['estimated_peptide_spot_count'].append(deconvolved_peptide.estimated_spot_count)
            data['hit_pool_ids'].append(';'.join([str(p) for p in deconvolved_peptide.hit_pool_ids]))
            data['hit_pools_count'].append(len(deconvolved_peptide.hit_pool_ids))
            data['deconvolution_result'].append(deconvolved_peptide.result.value)
            if self.pool_spot_counts is not None:
                hit_pool_spot_counts = []
                for hit_pool_id in deconvolved_peptide.hit_pool_ids:
                    hit_pool_spot_counts.append(self.pool_spot_counts[hit_pool_id])
                data['hit_pool_spot_counts'].append(';'.join([str(c) for c in hit_pool_spot_counts]))
            else:
                data['hit_pool_spot_counts'].append('')
            if self.plate_map is not None:
                if len(self.plate_map) > 0:
                    plate_well_ids = []
                    for hit_pool_id in deconvolved_peptide.hit_pool_ids:
                        plate_well = self.plate_map[hit_pool_id]
                        plate_well_ids.append('%i-%s' % (plate_well.plate_id, plate_well.well_id))
                    data['hit_plate_well_ids'].append(';'.join(plate_well_ids))
                else:
                    data['hit_plate_well_ids'].append('')
            else:
                data['hit_plate_well_ids'].append('')

        df = pd.DataFrame(data)
        df.sort_values(by=['peptide_id'], inplace=True)
        return df

    def __len__(self) -> int:
        return len(self.deconvolved_peptides)

    def __iter__(self):
        return iter(self.deconvolved_peptides)
