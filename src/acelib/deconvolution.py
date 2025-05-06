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


import numpy as np
import pandas as pd
from collections import defaultdict
from golfy import deconvolve
from typing import Dict, List
from .block_assignment import BlockAssignment
from .constants import DeconvolutionResult, DeconvolutionMethod
from .deconvolved_peptide import DeconvolvedPeptide
from .deconvolved_peptide_set import DeconvolvedPeptideSet
from .logger import get_logger
from .peptide import Peptide
from .utilities import convert_to_golfy_spot_counts


logger = get_logger(__name__)


def compute_background_spot_count(
        df_readout: pd.DataFrame,
        block_assignment: BlockAssignment,
        peptide_spot_counts: Dict[str,float],
        hit_peptide_ids: List[str]
) -> float:
    """
    Compute background spot count.

    Parameters:
        df_readout                  :   Pandas DataFrame with the following columns: 'pool_id', 'spot_count'.
        block_assignment            :   BlockAssignment object.
        peptide_spot_counts         :   Peptide spot counts (key = peptide ID, value = estimated spot).
        hit_peptide_ids             :   Hit peptide IDs (based on empirical deconvolution).

    Returns:
        background_spot_count       :   Background spot count.
    """
    deltas = []
    df_assignment = block_assignment.to_dataframe()
    for _, row in df_readout.iterrows():
        pool_id = row['pool_id']
        pool_spot_count = row['spot_count']
        hit_peptide_spot_counts_sum = 0.0
        non_hit_peptide_ids = []
        pool_peptide_ids = df_assignment.loc[df_assignment['pool_id'] == pool_id, 'peptide_id'].values.tolist()
        for peptide_id in pool_peptide_ids:
            if peptide_id in hit_peptide_ids:
                hit_peptide_spot_counts_sum += peptide_spot_counts[peptide_id]
            else:
                non_hit_peptide_ids.append(peptide_id)
        delta = pool_spot_count - hit_peptide_spot_counts_sum
        if len(non_hit_peptide_ids) == 0:
            delta = delta / len(pool_peptide_ids)
        else:
            delta = delta / len(non_hit_peptide_ids)
        deltas.append(delta)
    background_peptide_spot_count = float(np.mean(deltas))
    return background_peptide_spot_count


def perform_empirical_deconvolution(
        df_readout: pd.DataFrame,
        block_assignment: BlockAssignment,
        min_pool_spot_count: float,
        min_coverage: int
) -> DeconvolvedPeptideSet:
    """
    Perform empirical deconvolution of read-outs from an ELISpot experiment to identify hit peptides.

    Parameters:
        df_readout              :   Pandas DataFrame with the following columns: 'pool_id', 'spot_count'.
        block_assignment        :   BlockAssignment objects.
        min_pool_spot_count     :   Minimum pool spot count (inclusive).
        min_coverage            :   Minimum coverage.

    Returns:
        deconvolved_peptide_set :   DeconvolvedPeptideSet object.
    """
    # Step 1. Identify hit pool IDs
    df_assignment = block_assignment.to_dataframe()
    hit_pool_ids = df_readout.loc[df_readout['spot_count'] >= min_pool_spot_count,'pool_id'].unique()

    # Step 2. Identify hit pool IDs for each peptide ID
    # peptide_pools_dict:
    # {
    #     peptide_id_1: [pool_id_1,pool_id_2,...],
    #     ...
    # }
    peptide_hit_pools_dict = defaultdict(list)
    for pool_id in hit_pool_ids:
        hit_peptide_ids = df_assignment.loc[df_assignment['pool_id'] == pool_id, 'peptide_id'].values.tolist()
        for peptide_id in hit_peptide_ids:
            peptide_hit_pools_dict[peptide_id].append(pool_id)
    for peptide_id in df_assignment['peptide_id'].unique():
        if peptide_id not in peptide_hit_pools_dict.keys():
            peptide_hit_pools_dict[peptide_id] = []

    # Step 3. Identify hit peptides (peptides that meet the minimum coverage)
    hit_peptide_ids = set()
    for peptide_id, pool_ids in peptide_hit_pools_dict.items():
        if len(pool_ids) >= min_coverage:
            hit_peptide_ids.add(peptide_id)

    # Step 4. Identify hit pool IDs and the associated peptide IDs
    # hit_pool_peptides_dict:
    # {
    #     pool_id_1: [peptide_id_1,peptide_id_2,...],
    #     ...
    # }
    hit_pool_peptides_dict = defaultdict(list)
    for hit_peptide_id in hit_peptide_ids:
        pool_ids = peptide_hit_pools_dict[hit_peptide_id]
        for pool_id in pool_ids:
            hit_pool_peptides_dict[int(pool_id)].append(hit_peptide_id)

    # Step 5. Identify candidate peptides (second-round assay peptides) for each hit peptide
    candidate_peptide_ids = set()
    for hit_peptide_id in hit_peptide_ids:
        pool_ids = df_assignment.loc[df_assignment['peptide_id'] == hit_peptide_id, 'pool_id'].values.tolist()
        unique = False
        for pool_id in pool_ids:
            if len(hit_pool_peptides_dict[int(pool_id)]) == 1:
                unique = True
        if not unique:
            candidate_peptide_ids.add(hit_peptide_id)

    # Step 6. Prepare output
    pool_spot_counts = {}
    for pool_id in df_readout['pool_id'].unique():
        pool_spot_count = df_readout.loc[df_readout['pool_id'] == pool_id,'spot_count'].values.tolist()[0]
        pool_spot_counts[pool_id] = pool_spot_count

    deconvolved_peptide_set = DeconvolvedPeptideSet(
        deconvolution_method=DeconvolutionMethod.EMPIRICAL,
        pool_spot_counts=pool_spot_counts,
        min_pool_spot_count=min_pool_spot_count,
        min_coverage=min_coverage
    )

    for peptide_id in df_assignment['peptide_id'].unique():
        peptide_sequence = df_assignment.loc[df_assignment['peptide_id'] == peptide_id ,'peptide_sequence'].values[0]
        peptide = Peptide(id=peptide_id, sequence=peptide_sequence)
        hit_pool_ids = peptide_hit_pools_dict[peptide_id]
        estimated_spot_count = len(hit_pool_ids)
        if peptide_id in hit_peptide_ids:
            if peptide_id in candidate_peptide_ids:
                deconvolved_peptide = DeconvolvedPeptide(
                    peptide=peptide,
                    estimated_spot_count=estimated_spot_count,
                    result=DeconvolutionResult.CANDIDATE_HIT,
                    hit_pool_ids=hit_pool_ids
                )
                deconvolved_peptide_set.add(deconvolved_peptide=deconvolved_peptide)
            else:
                deconvolved_peptide = DeconvolvedPeptide(
                    peptide=peptide,
                    estimated_spot_count=estimated_spot_count,
                    result=DeconvolutionResult.CONFIDENT_HIT,
                    hit_pool_ids=hit_pool_ids
                )
                deconvolved_peptide_set.add(deconvolved_peptide=deconvolved_peptide)
        else:
            deconvolved_peptide = DeconvolvedPeptide(
                peptide=peptide,
                estimated_spot_count=estimated_spot_count,
                result=DeconvolutionResult.NOT_A_HIT,
                hit_pool_ids=hit_pool_ids
            )
            deconvolved_peptide_set.add(deconvolved_peptide=deconvolved_peptide)

    deconvolved_peptide_set.plate_map = block_assignment.plate_map

    return deconvolved_peptide_set


def perform_statistical_deconvolution(
        df_readout: pd.DataFrame,
        block_assignment: BlockAssignment,
        method: DeconvolutionMethod,
        min_peptide_spot_count: float,
        verbose: bool = True
) -> DeconvolvedPeptideSet:
    """
    Perform statistical deconvolution.

    Parameters:
        df_readout              :   Pandas DataFrame with the following columns: 'pool_id', 'spot_count'.
        block_assignment        :   BlockAssignment objects.
        method                  :   Deconvolution method.
        min_peptide_spot_count  :   Minimum peptide spot count.
        verbose                 :   If True, print logs.

    Returns:
        deconvolved_peptide_set :   DeconvolvedPeptideSet object.
    """
    assert method not in [DeconvolutionMethod.EMPIRICAL, DeconvolutionMethod.CONSTRAINED_EM]

    # Step 1. Run golfy
    golfy_design, peptide_indices = block_assignment.to_golfy_design()
    spot_counts = {}
    for _, row in df_readout.iterrows():
        spot_counts[int(row['pool_id'])] = int(row['spot_count'])
    golfy_spot_counts = convert_to_golfy_spot_counts(
        spot_counts=spot_counts,
        block_assignment=block_assignment
    )
    golfy_deconvolve_result = deconvolve(
        s=golfy_design,
        spot_counts=golfy_spot_counts,
        method=method.value,
        min_peptide_activity=min_peptide_spot_count,
        verbose=verbose
    )

    # Step 2. Prepare output
    pool_spot_counts = {}
    for pool_id in df_readout['pool_id'].unique():
        pool_spot_count = df_readout.loc[df_readout['pool_id'] == pool_id,'spot_count'].values.tolist()[0]
        pool_spot_counts[pool_id] = pool_spot_count

    deconvolved_peptide_set = DeconvolvedPeptideSet(
        deconvolution_method=method,
        min_peptide_spot_count=min_peptide_spot_count,
        pool_spot_counts=pool_spot_counts
    )

    for peptide_idx,peptide_activity in enumerate(golfy_deconvolve_result.activity_per_peptide):
        peptide_id = peptide_indices[peptide_idx]
        peptide_sequence = block_assignment.get_peptide_sequence(peptide_id=peptide_id)
        if peptide_idx in golfy_deconvolve_result.high_confidence_hits:
            deconvolution_result = DeconvolutionResult.CONFIDENT_HIT
        else:
            deconvolution_result = DeconvolutionResult.NOT_A_HIT
        peptide = Peptide(id=peptide_id, sequence=peptide_sequence)
        deconvolved_peptide = DeconvolvedPeptide(
            peptide=peptide,
            estimated_spot_count=peptide_activity,
            result=deconvolution_result,
            hit_pool_ids=[]
        )
        deconvolved_peptide_set.add(deconvolved_peptide=deconvolved_peptide)

    deconvolved_peptide_set.plate_map = block_assignment.plate_map

    return deconvolved_peptide_set
