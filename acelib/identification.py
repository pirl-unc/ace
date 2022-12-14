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
The purpose of this python3 script is to
implement functions related to identification of positive peptide hits.
"""


from __future__ import print_function, division, absolute_import


import pandas as pd
from .logger import get_logger
from .default_parameters import *


logger = get_logger(__name__)


def identify_hit_peptides(df_readout: pd.DataFrame,
                          df_configuration: pd.DataFrame,
                          min_positive_spot_count: int = MIN_POSITIVE_SPOT_COUNT) -> pd.DataFrame:
    """
    Identifies hit peptide IDs given read-outs from an ELIspot experiment.

    Parameters
    ----------
    df_readout              :   DataFrame of ELIspot experiment read-out results.
                                Expected columns:
                                'pool_id'
                                'spot_count'
    df_configuration        :   DataFrame of ELIspot configuration.
                                Expected columns:
                                'coverage_id'
                                'pool_id'
                                'peptide_id'
    min_positive_spot_count :   Minimum number of spots for a pool to be
                                considered a positive hit.

    Returns
    -------
    df_hits                 :   DataFrame with the following columns:
                                'hit_peptide_id'
                                'num_coverage'
                                'pool_ids'
    """
    # Step 1. Figure out the coverage
    max_coverage = max([int(i.split('_')[1]) for i in df_configuration['coverage_id'].unique()])

    # Step 2. Identify hit pool IDs
    hit_pool_ids = df_readout.loc[df_readout['spot_count'] >= min_positive_spot_count, 'pool_id'].unique()

    # Step 3. Identify hit peptide IDs
    hit_peptide_ids = set()
    for curr_pool_id in hit_pool_ids:
        for curr_peptide_id in df_configuration.loc[df_configuration['pool_id'] == curr_pool_id, 'peptide_id'].values.tolist():
            hit_peptide_ids.add(curr_peptide_id)

    # Step 4. Identify coverage and pool IDs for each hit peptide
    data = {
        'hit_peptide_id': [],
        'num_coverage': [],
        'pool_ids': []
    }
    for curr_hit_peptide_id in hit_peptide_ids:
        curr_hit_peptide_pool_ids = []
        for curr_pool_id in hit_pool_ids:
            curr_pool_peptide_ids = df_configuration.loc[df_configuration['pool_id'] == curr_pool_id, 'peptide_id'].values.tolist()
            if curr_hit_peptide_id in curr_pool_peptide_ids:
                curr_hit_peptide_pool_ids.append(curr_pool_id)
        data['hit_peptide_id'].append(curr_hit_peptide_id)
        data['num_coverage'].append(len(curr_hit_peptide_pool_ids))
        data['pool_ids'].append(','.join(curr_hit_peptide_pool_ids))

    df_hits = pd.DataFrame(data)


    #while len(df_hits[df_hits['num_coverage'] == max_coverage]) == 0:
    #    max_coverage -=1

    #df_hits = df_hits[df_hits['num_coverage'] == max_coverage]

    if len(df_hits) == 0:
        logger.warning("No peptide hit was identified.")

    return df_hits


