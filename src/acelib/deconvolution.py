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


import golfy
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List
from .block_assignment import BlockAssignment
from .constants import DeconvolutionLabels
from .logger import get_logger
from .types import *


logger = get_logger(__name__)


@dataclass(frozen=True)
class DeconvolutionResult:
    hit_peptides: HitPeptides = field(default_factory=list)

    def add_peptide(
            self,
            peptide_id: str,
            peptide_sequence: str,
            estimated_peptide_spot_count: float,
            label: str,
            hit_pool_ids: List[PoolId]
    ):
        """
        Add a peptide.

        Parameters
        ----------
        peptide_id                      :   Peptide ID.
        peptide_sequence                :   Peptide sequence.
        estimated_peptide_spot_count    :   Estimated peptide spot count.
        label                           :   Deconvolution label.
        hit_pool_ids                    :   Hit pool IDs.
        """
        self.hit_peptides.append(
            (peptide_id, 
             peptide_sequence, 
             estimated_peptide_spot_count,
             label,
             hit_pool_ids)
        )

    def to_dataframe(self):
        data = {
            'peptide_id': [],
            'peptide_sequence': [],
            'estimated_peptide_spot_count': [],
            'hit_pool_ids': [],
            'hit_pools_count': [],
            'deconvolution_result': []
        }
        for peptide_id, peptide_sequence, peptide_spot_count, label, pool_ids in self.hit_peptides:
            data['peptide_id'].append(peptide_id)
            data['peptide_sequence'].append(peptide_sequence)
            data['estimated_peptide_spot_count'].append(peptide_spot_count)
            data['hit_pool_ids'].append(';'.join([str(p) for p in pool_ids]))
            data['hit_pools_count'].append(len(pool_ids))
            data['deconvolution_result'].append(label)
        df = pd.DataFrame(data)
        df.sort_values(by=['peptide_id'], inplace=True)
        return df


def convert_to_golfy_spotcounts(
        spot_counts: SpotCounts, 
        block_assignment: BlockAssignment
) -> golfy.SpotCounts:
    """
    Converts spot counts to golfy SpotCounts.

    Parameters
    ----------
    spot_counts             :   Dictionary where
                                key     = pool ID
                                value   = spot count
    block_assignment        :   BlockAssignment object.

    Returns
    -------
    spot_counts_            :   golfy.SpotCounts.
    """
    spot_counts_ = {}
    for coverage in block_assignment.assignments.keys():
        spot_counts_[coverage-1] = {}
        for pool in block_assignment.assignments[coverage].keys():
            spot_counts_[coverage-1][pool-1] = spot_counts[pool]
    return spot_counts_


def empirical_deconvolve(
        hit_pool_ids: List[PoolId],
        df_assignment: pd.DataFrame,
        min_coverage: int
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
    min_coverage            :   Minimum coverage.

    Returns
    -------
    deconvolution_result    :   EmpiricalDeconvolutionResult object.
    """
    # Step 1. Identify hit peptide IDs
    # key   = peptide ID
    # value = pool IDs
    hit_peptides_dict = defaultdict(list)
    for curr_pool_id in hit_pool_ids:
        curr_hit_peptide_ids = df_assignment.loc[df_assignment['pool_id'] == curr_pool_id, 'peptide_id'].values.tolist()
        for curr_hit_peptide_id in curr_hit_peptide_ids:
            hit_peptides_dict[curr_hit_peptide_id].append(curr_pool_id)
    for peptide_id in df_assignment['peptide_id'].unique():
        if peptide_id not in hit_peptides_dict.keys():
            hit_peptides_dict[peptide_id] = []

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

    # Step 3. Identify peptide coverage
    df_hits_ = df_hits[df_hits['num_coverage'] >= min_coverage]
    if len(df_hits_) == 0:
        logger.info('Returning as there are no peptides with the desired minimum hit coverage (%ix).' % min_coverage)
        deconvolution_result = DeconvolutionResult()
        for index, row in df_hits.iterrows():
            peptide_id = row['peptide_id']
            peptide_sequence = df_assignment.loc[df_assignment['peptide_id'] == peptide_id, 'peptide_sequence'].values[0]
            if row['pool_ids'] != '':
                pool_ids = [int(pool_id) for pool_id in row['pool_ids'].split(';')]
            else:
                pool_ids = []
            deconvolution_result.add_peptide(
                peptide_id=peptide_id,
                peptide_sequence=peptide_sequence,
                estimated_peptide_spot_count=row['num_coverage'],
                label=DeconvolutionLabels.NOT_A_HIT,
                hit_pool_ids=pool_ids
            )
        return deconvolution_result

    # Step 4. Identify hit pool IDs and the associated peptide IDs
    hit_pool_ids_dict = defaultdict(list)  # key = pool ID, value = list of peptide IDs
    for _, value in df_hits_.iterrows():
        peptide_id = value['peptide_id']
        pool_ids = value['pool_ids'].split(';')
        for pool_id in pool_ids:
            hit_pool_ids_dict[int(pool_id)].append(peptide_id)

    # Step 5. For the peptides that have the maximum coverage, identify second-round assay peptides
    second_round_assay_peptide_ids = set()
    for peptide_id in df_hits_['peptide_id'].unique():
        pool_ids = df_assignment.loc[df_assignment['peptide_id'] == peptide_id, 'pool_id'].values.tolist()
        unique = False
        for pool_id in pool_ids:
            if len(hit_pool_ids_dict[int(pool_id)]) == 1:
                unique = True
        if not unique:
            second_round_assay_peptide_ids.add(peptide_id)

    # Step 6. Deconvolve hit peptide IDs
    data = {
        'peptide_id': [],
        'peptide_sequence': [],
        'deconvolution_label': []
    }
    for peptide_id in df_assignment['peptide_id'].unique():
        peptide_sequence = df_assignment.loc[df_assignment['peptide_id'] == peptide_id ,'peptide_sequence'].values[0]
        if peptide_id in df_hits_['peptide_id'].unique():
            if peptide_id not in second_round_assay_peptide_ids:
                data['peptide_id'].append(peptide_id)
                data['peptide_sequence'].append(peptide_sequence)
                data['deconvolution_label'].append(DeconvolutionLabels.CONFIDENT_HIT)
            else:
                data['peptide_id'].append(peptide_id)
                data['peptide_sequence'].append(peptide_sequence)
                data['deconvolution_label'].append(DeconvolutionLabels.CANDIDATE_HIT)
        else:
            data['peptide_id'].append(peptide_id)
            data['peptide_sequence'].append(peptide_sequence)
            data['deconvolution_label'].append(DeconvolutionLabels.NOT_A_HIT)
    df_deconvolution = pd.DataFrame(data)
    df_deconvolution = pd.merge(df_hits, df_deconvolution, on=['peptide_id'])
    
    assert len(df_assignment['peptide_id'].unique()) == len(df_deconvolution), 'All peptide IDs must be present in the deconvolution result.'

    deconvolution_result = DeconvolutionResult()
    for _, row in df_deconvolution.iterrows():
        if row['pool_ids'] == '':
            pool_ids = []
        else:
            pool_ids = [int(p) for p in row['pool_ids'].split(';')]
        deconvolution_result.add_peptide(
            peptide_id=row['peptide_id'],
            peptide_sequence=row['peptide_sequence'],
            estimated_peptide_spot_count=len(pool_ids),
            label=row['deconvolution_label'],
            hit_pool_ids=pool_ids
        )
    return deconvolution_result

