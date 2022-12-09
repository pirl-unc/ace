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
implement functions related to verifying an ELIspot configuration.
"""


import pandas as pd
from .logger import get_logger


logger = get_logger(__name__)


def verify_configuration_constraints(df_configuration: pd.DataFrame,
                                     n_peptides: int,
                                     n_peptides_per_pool: int,
                                     n_coverage: int) -> bool:
    """
    Verifies whether a given ELIspot configuration satisfies the following constraints:
    1. Total number of peptides is 'n_peptides'.
    2. Each pool contains 'n_peptides_per_pool' number of peptides.
    3. Each coverage has (n_peptides / n_peptides_per_pool) number of pools.
    4. Each peptide belong in exactly one unique combination of pool IDs.

    Parameters
    ----------
    df_configuration        :   DataFrame of ELIspot configuration.
                                Expected columns:
                                'coverage_id'
                                'pool_id'
                                'peptide_id'
    n_peptides              :   Total number of peptides.
    n_peptides_per_pool     :   Number of peptides per pool.
    n_coverage              :   Coverage.

    Returns
    -------
    True if the input configuration meets all desired constraints.
    False otherwise.
    """
    # Step 1. Check that the total number of peptides is 'n_peptides'
    n_peptides_configuration = len(df_configuration['peptide_id'].unique())
    if n_peptides_configuration != n_peptides:
        logger.info("%i total unique peptides in the configuration. "
                    "Input total: %i. "
                    "The total number of peptides does not match the input number." %
                    (n_peptides_configuration, n_peptides))
        return False

    # Step 2. Check that each pool contains 'n_peptides_per_pool' number of peptides
    for curr_pool_id in df_configuration['pool_id'].unique():
        curr_pool_peptide_ids = df_configuration.loc[
            df_configuration['pool_id'] == curr_pool_id,
            'peptide_id'
        ].unique()
        if len(curr_pool_peptide_ids) != n_peptides_per_pool:
            logger.info("%i unique peptides in %s. "
                        "There are supposed to be %i." %
                        (len(curr_pool_peptide_ids),
                         curr_pool_id,
                         n_peptides_per_pool))
            return False

    # Step 3. Check that in each coverage, there are (n_peptides / n_peptides_per_pool) number of pools
    n_pools_per_coverage = int(n_peptides / n_peptides_per_pool)
    for curr_coverage_id in df_configuration['coverage_id'].unique():
        curr_coverage_pool_ids = df_configuration.loc[
            df_configuration['coverage_id'] == curr_coverage_id,
            'pool_id'
        ].unique()
        if len(curr_coverage_pool_ids) != n_pools_per_coverage:
            logger.info("%i unique pools in %s. "
                        "There are supposed to be %i." %
                        (len(curr_coverage_pool_ids),
                         curr_coverage_id,
                         n_pools_per_coverage))
            return False

    # Step 4. Check that each peptide belongs in exactly one unique combination of pool IDs
    pool_id_combinations = set()
    for curr_peptide_id in df_configuration['peptide_id'].unique():
        curr_peptide_pool_ids = df_configuration.loc[
            df_configuration['peptide_id'] == curr_peptide_id,
            'pool_id'
        ].unique()
        if len(curr_peptide_pool_ids) != n_coverage:
            logger.info("%s assigned in %i unique pools. "
                        "Coverage is supposed to be %i." %
                        (curr_peptide_id,
                         len(curr_peptide_pool_ids),
                         n_coverage))
            return False
        curr_peptide_pool_ids = sorted(curr_peptide_pool_ids)
        pool_id_combinations.add(','.join(curr_peptide_pool_ids))
    if len(pool_id_combinations) != n_peptides:
        logger.info("%i unique combinations of pool IDs. "
                    "Total number of peptides: %i. "
                    "Each peptide is supposed to be assigned a unique combination of pools." %
                    (len(pool_id_combinations), n_peptides))
        return False

    # Meets all constraints
    logger.info("Configuration meets all ACE constraints.")
    return True
