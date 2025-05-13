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
import numpy as np
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
    min_pool_spot_count: float
    min_coverage: int
    deconvolution_method: DeconvolutionMethod
    pool_spot_counts: Dict[int,float] # key = pool ID, value = pool spot count
    min_peptide_spot_count: Optional[float] = None
    background_spot_count: Optional[float] = None
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

    def get_pool_spot_count(self, pool_id: int) -> float:
        return self.pool_spot_counts[pool_id]

    def get_plate_well_spot_count(self, plate_id: int, well_id: str) -> float:
        for pool_id, plate_well in self.plate_map.items():
            if plate_well.plate_id == plate_id and plate_well.well_id == well_id:
                return self.get_pool_spot_count(pool_id=pool_id)
        raise Exception('Could not find the spot count for %i-%s' % (plate_id, well_id))

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
                                - 'hit_plate_well_ids'
                                - 'hit_plate_well_ids_count'
                                - 'hit_plate_well_spot_counts'
                                - 'hit_plate_well_spot_counts_average'
                                - 'hit_plate_well_spot_counts_standard_deviation'
                                - 'hit_plate_well_spot_counts_variation_coefficient'
                                - 'hit_pool_ids'
                                - 'hit_pool_ids_count'
                                - 'hit_pool_spot_counts'
                                - 'hit_pool_spot_counts_average'
                                - 'hit_pool_spot_counts_standard_deviation'
                                - 'hit_pool_spot_counts_variation_coefficient'
                                - 'estimated_peptide_spot_count'
                                - 'deconvolution_result'
                                - 'hit_plate_well_spot_count_1'
                                - ...
                                - 'hit_pool_spot_count_1'
                                - ...
        """
        data = {
            'peptide_id': [],
            'peptide_sequence': [],
            'hit_plate_well_ids': [],
            'hit_plate_well_ids_count': [],
            'hit_plate_well_spot_counts': [],
            'hit_plate_well_spot_counts_average': [],
            'hit_plate_well_spot_counts_standard_deviation': [],
            'hit_plate_well_spot_counts_variation_coefficient': [],
            'hit_pool_ids': [],
            'hit_pool_ids_count': [],
            'hit_pool_spot_counts': [],
            'hit_pool_spot_counts_average': [],
            'hit_pool_spot_counts_standard_deviation': [],
            'hit_pool_spot_counts_variation_coefficient': [],
            'estimated_peptide_spot_count': [],
            'deconvolution_result': []
        }

        for deconvolved_peptide in self.deconvolved_peptides:
            data['peptide_id'].append(deconvolved_peptide.id)
            data['peptide_sequence'].append(deconvolved_peptide.sequence)

            if self.plate_map:
                plate_well_ids = []
                plate_well_spot_counts = []
                for hit_pool_id in deconvolved_peptide.hit_pool_ids:
                    plate_well = self.plate_map[hit_pool_id]
                    spot_count = self.get_plate_well_spot_count(plate_id=plate_well.plate_id, well_id=plate_well.well_id)
                    plate_well_ids.append('%i-%s' % (plate_well.plate_id, plate_well.well_id))
                    plate_well_spot_counts.append(spot_count)
                data['hit_plate_well_ids'].append(';'.join(plate_well_ids))
                data['hit_plate_well_ids_count'].append(len(plate_well_ids))
                data['hit_plate_well_spot_counts'].append(';'.join(map(str, plate_well_spot_counts)))

                for i in range(0,self.min_coverage):
                    if i < len(plate_well_spot_counts):
                        spot_count = plate_well_spot_counts[i]
                    else:
                        spot_count = ''
                    if 'hit_plate_well_spot_count_%i' % (i+1) not in data.keys():
                        data['hit_plate_well_spot_count_%i' % (i+1)] = [spot_count]
                    else:
                        data['hit_plate_well_spot_count_%i' % (i+1)].append(spot_count)

                if len(plate_well_spot_counts) > 0:
                    counts_array = np.array(plate_well_spot_counts, dtype=float)
                    spot_count_average = float(np.mean(counts_array))
                    spot_count_std = float(np.std(counts_array))
                    if spot_count_average != 0:
                        spot_count_var_coeff = spot_count_std / spot_count_average
                    else:
                        spot_count_var_coeff = ''
                else:
                    spot_count_average = ''
                    spot_count_std = ''
                    spot_count_var_coeff = ''
                data['hit_plate_well_spot_counts_average'].append(spot_count_average)
                data['hit_plate_well_spot_counts_standard_deviation'].append(spot_count_std)
                data['hit_plate_well_spot_counts_variation_coefficient'].append(spot_count_var_coeff)
            else:
                for key in [
                    'hit_plate_well_ids',
                    'hit_plate_well_ids_count',
                    'hit_plate_well_spot_counts',
                    'hit_plate_well_spot_counts_average',
                    'hit_plate_well_spot_counts_standard_deviation',
                    'hit_plate_well_spot_counts_variation_coefficient'
                ]:
                    data[key].append('')

            if self.pool_spot_counts:
                pool_ids = []
                pool_spot_counts = []
                for hit_pool_id in deconvolved_peptide.hit_pool_ids:
                    spot_count = self.get_pool_spot_count(pool_id=hit_pool_id)
                    pool_ids.append(hit_pool_id)
                    pool_spot_counts.append(spot_count)
                data['hit_pool_ids'].append(';'.join(map(str, pool_ids)))
                data['hit_pool_ids_count'].append(len(pool_ids))
                data['hit_pool_spot_counts'].append(';'.join(map(str, pool_spot_counts)))

                for i in range(0,self.min_coverage):
                    if i < len(pool_spot_counts):
                        spot_count = pool_spot_counts[i]
                    else:
                        spot_count = ''
                    if 'hit_pool_spot_count_%i' % (i+1) not in data.keys():
                        data['hit_pool_spot_count_%i' % (i+1)] = [spot_count]
                    else:
                        data['hit_pool_spot_count_%i' % (i+1)].append(spot_count)

                if len(pool_spot_counts) > 0:
                    counts_array = np.array(pool_spot_counts, dtype=float)
                    spot_count_average = float(np.mean(counts_array))
                    spot_count_std = float(np.std(counts_array))
                    if spot_count_average != 0:
                        spot_count_var_coeff = spot_count_std / spot_count_average
                    else:
                        spot_count_var_coeff = ''
                else:
                    spot_count_average = ''
                    spot_count_std = ''
                    spot_count_var_coeff = ''
                data['hit_pool_spot_counts_average'].append(spot_count_average)
                data['hit_pool_spot_counts_standard_deviation'].append(spot_count_std)
                data['hit_pool_spot_counts_variation_coefficient'].append(spot_count_var_coeff)
            else:
                for key in [
                    'hit_pool_ids',
                    'hit_pool_ids_count',
                    'hit_pool_spot_counts',
                    'hit_pool_spot_counts_average',
                    'hit_pool_spot_counts_standard_deviation',
                    'hit_pool_spot_counts_variation_coefficient'
                ]:
                    data[key].append('')

            data['estimated_peptide_spot_count'].append(deconvolved_peptide.estimated_spot_count)
            data['deconvolution_result'].append(deconvolved_peptide.result.value)

        df = pd.DataFrame(data)
        df.sort_values(by=['peptide_id'], inplace=True)
        return df

    def __len__(self) -> int:
        return len(self.deconvolved_peptides)

    def __iter__(self):
        return iter(self.deconvolved_peptides)
