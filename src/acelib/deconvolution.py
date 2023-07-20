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
The purpose of this python3 script is to implement functions related to deconvolution.
"""


from collections import defaultdict
from dataclasses import dataclass, field
from typing import List
from .block_assignment import BlockAssignment
from .block_design import BlockDesign
from .constants import DeconvolutionLabels
from .logger import get_logger
from .types import *


logger = get_logger(__name__)


@dataclass(frozen=True)
class DeconvolutionResult:
    hit_peptides: HitPeptides = field(default_factory=list)

    def add_hit_peptide(
            self,
            peptide_id: str,
            peptide_sequence: str,
            deconvolution_label: str,
            pool_ids: List[PoolId]
    ):
        """
        Add a hit peptide.

        Parameters
        ----------
        peptide_id              :   Peptide ID.
        peptide_sequence        :   Peptide sequence.
        deconvolution_label     :   Deconvolution label.
        pool_ids                :   Pool IDs.
        """
        self.hit_peptides.append((peptide_id, peptide_sequence, deconvolution_label, pool_ids))

    def to_dataframe(self):
        data = {
            'peptide_id': [],
            'peptide_sequence': [],
            'pool_ids': [],
            'num_coverage': [],
            'deconvolution_result': []
        }
        for peptide_id, peptide_sequence, label, pool_ids in self.hit_peptides:
            data['peptide_id'].append(peptide_id)
            data['peptide_sequence'].append(peptide_sequence)
            data['pool_ids'].append(';'.join([str(p) for p in pool_ids]))
            data['num_coverage'].append(len(pool_ids))
            data['deconvolution_result'].append(label)
        return pd.DataFrame(data)


def deconvolve_hit_peptides(
        hit_pool_ids: List[PoolId],
        df_assignment: pd.DataFrame,
        num_coverage: int
) -> DeconvolutionResult:
    """
    Identifies hit peptide IDs given read-outs from an ELISpot experiment.

    Parameters
    ----------
    hit_pool_ids            :   Hit pool IDs.
    df_assignment           :   pd.DataFrame with the following columns:
                                'coverage_id'
                                'pool_id'
                                'peptide_id'
                                'peptide_sequence'
    num_coverage            :   Coverage.

    Returns
    -------
    deconvolution_result    :   DeconvolutionResult object.
    """
    # Step 1. Deconvolve hit peptide IDs
    # key = peptide ID
    # value = pool IDs
    hit_peptides_dict = defaultdict(list)
    for curr_pool_id in hit_pool_ids:
        curr_hit_peptide_ids = df_assignment.loc[df_assignment['pool_id'] == curr_pool_id, 'peptide_id'].values.tolist()
        for curr_hit_peptide_id in curr_hit_peptide_ids:
            hit_peptides_dict[curr_hit_peptide_id].append(curr_pool_id)

    # Step 2. Identify coverage and pool IDs for each hit peptide
    data = {
        'peptide_id': [],
        'pool_ids': [],
        'num_coverage': []
    }
    for key, value in hit_peptides_dict.items():
        data['peptide_id'].append(key)
        data['pool_ids'].append(';'.join([str(i) for i in value]))
        data['num_coverage'].append(len(value))
    df_hits = pd.DataFrame(data)

    # Step 3. Identify peptide maximum coverage
    df_hits = df_hits[df_hits['num_coverage'] == num_coverage]
    if len(df_hits) == 0:
        logger.info('Returning as there are no peptides with the desired hit coverage (%ix).' % num_coverage)
        return DeconvolutionResult()

    # Step 4. Identify hit pool IDs and the associated peptide IDs
    hit_pool_ids_dict = defaultdict(list)  # key = pool ID, value = list of peptide IDs
    for index, value in df_hits.iterrows():
        peptide_id = value['peptide_id']
        pool_ids = value['pool_ids'].split(';')
        for pool_id in pool_ids:
            hit_pool_ids_dict[int(pool_id)].append(peptide_id)

    # Step 5. For the peptides that have the maximum coverage,
    # identify second-round assay peptides
    second_round_assay_peptide_ids = set()
    for peptide_id in df_hits['peptide_id'].unique():
        pool_ids = df_assignment.loc[df_assignment['peptide_id'] == peptide_id, 'pool_id'].values.tolist()
        unique = False
        for pool_id in pool_ids:
            if len(hit_pool_ids_dict[int(pool_id)]) == 1:
                unique = True
        if not unique:
            second_round_assay_peptide_ids.add(peptide_id)

    # Step 6. Deconvolve hit peptide IDs.
    data = {
        'peptide_id': [],
        'peptide_sequence': [],
        'deconvolution_label': []
    }
    for peptide_id in df_hits['peptide_id'].unique():
        peptide_sequence = df_assignment.loc[df_assignment['peptide_id'] == peptide_id ,'peptide_sequence'].values[0]
        if peptide_id not in second_round_assay_peptide_ids:
            data['peptide_id'].append(peptide_id)
            data['peptide_sequence'].append(peptide_sequence)
            data['deconvolution_label'].append(DeconvolutionLabels.HIT)
        else:
            data['peptide_id'].append(peptide_id)
            data['peptide_sequence'].append(peptide_sequence)
            data['deconvolution_label'].append(DeconvolutionLabels.CANDIDATE_HIT)
    df_deconvolution = pd.DataFrame(data)
    df_deconvolution = pd.merge(df_hits, df_deconvolution, on=['peptide_id'])

    deconvolution_result = DeconvolutionResult()
    for index, row in df_deconvolution.iterrows():
        deconvolution_result.add_hit_peptide(
            peptide_id=row['peptide_id'],
            peptide_sequence=row['peptide_sequence'],
            deconvolution_label=row['deconvolution_label'],
            pool_ids=[int(p) for p in row['pool_ids'].split(';')]
        )

    return deconvolution_result

