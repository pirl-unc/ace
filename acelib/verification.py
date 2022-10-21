#!/usr/bin/python3

"""
The purpose of this python3 script is to implement functions related to
verifying configuration and identifying hit peptide IDs.

Author: Jin Seok (Andy) Lee, Dhuvarakesh Karthikeyan

Last updated date: Aug 15, 2022
"""


import pandas as pd
from acelib.logging import get_logger


logger = get_logger(__name__)


def identify_hit_peptides(df_configuration, hit_pool_ids):
    """
    Identifies hit peptides based on hit pool IDs.

    Args
    ----
    df_configuration    :   DataFrame with the following columns:
                            'coverage_id', 'pool_id', 'peptide_id'
    hit_pool_ids        :   List of pool IDs.

    Returns
    -------
    DataFrame with the following columns:
    'hit_peptide_id', 'num_coverage', 'pool_ids'
    """
    # Step 1. Figure out the coverage
    n_coverage = max([int(i.split('_')[1]) for i in df_configuration['coverage_id'].unique()])

    # Step 2. Identify hit peptide IDs
    hit_peptide_ids = set()
    for curr_pool_id in hit_pool_ids:
        for curr_peptide_id in df_configuration.loc[df_configuration['pool_id'] == curr_pool_id, 'peptide_id'].values.tolist():
            hit_peptide_ids.add(curr_peptide_id)

    # Step 3. Identify coverage and pool IDs for each hit peptide
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

    df = pd.DataFrame(data)
    return df


def verify_configuration_constraints(df_configuration, n_peptides, n_peptides_per_pool, n_coverage):
    """
    Verifies the following constraints for a configuration:
    1. Total number of peptides is 'n_peptides'.
    2. Each pool contains 'n_peptides_per_pool' number of peptides.
    3. Each coverage has (n_peptides / n_peptides_per_pool) number of pools.
    4. Each peptide belong in exactly one unique combination of pool IDs.

    Args
    ----
    df_configuration        :   DataFrame with the following columns:
                                'coverage_id', 'pool_id', 'peptide_id'
    n_peptides              :   Total number of peptides.
    n_peptides_per_pool     :   Number of peptides per pool.
    n_coverage              :   Coverage.

    Returns
    -------
    True if the input configuration meets all desired constraints
    for a bona fide configuration. False otherwise.
    """
    # Step 1. Check that the total number of peptides is 'n_peptides'
    n_peptides_configuration = len(df_configuration['peptide_id'].unique())
    if n_peptides_configuration != n_peptides:
        logger.info("%i total unique peptides in the configuration. Input total: %i. The total number of peptides does not match the input number." % (n_peptides_configuration, n_peptides))
        return False

    # Step 2. Check that each pool contains 'n_peptides_per_pool' number of peptides
    for curr_pool_id in df_configuration['pool_id'].unique():
        curr_pool_peptide_ids = df_configuration.loc[df_configuration['pool_id'] == curr_pool_id, 'peptide_id'].unique()
        if len(curr_pool_peptide_ids) != n_peptides_per_pool:
            logger.info("%i unique peptides in %s. There are supposed to be %i." % (len(curr_pool_peptide_ids), curr_pool_id, n_peptides_per_pool))
            return False

    # Step 3. Check that in each coverage, there are (n_peptides / n_peptides_per_pool) number of pools
    n_pools_per_coverage = int(n_peptides / n_peptides_per_pool)
    for curr_coverage_id in df_configuration['coverage_id'].unique():
        curr_coverage_pool_ids = df_configuration.loc[df_configuration['coverage_id'] == curr_coverage_id, 'pool_id'].unique()
        if len(curr_coverage_pool_ids) != n_pools_per_coverage:
            logger.info("%i unique pools in %s. There are supposed to be %i." % (len(curr_coverage_pool_ids), curr_coverage_id, n_pools_per_coverage))
            return False

    # Step 4. Check that each peptide belongs in exactly one unique combination of pool IDs
    pool_id_combinations = set()
    for curr_peptide_id in df_configuration['peptide_id'].unique():
        curr_peptide_pool_ids = df_configuration.loc[df_configuration['peptide_id'] == curr_peptide_id, 'pool_id'].unique()
        if len(curr_peptide_pool_ids) != n_coverage:
            logger.info("%s assigned in %i unique pools. Coverage is supposed to be %i." % (curr_peptide_id, len(curr_peptide_pool_ids), n_coverage))
            return False
        curr_peptide_pool_ids = sorted(curr_peptide_pool_ids)
        pool_id_combinations.add(','.join(curr_peptide_pool_ids))
    if len(pool_id_combinations) != n_peptides:
        logger.info("%i unique combinations of pool IDs. Total number of peptides: %i. Each peptide is supposed to be assigned a unique combination of pools." % (len(pool_id_combinations), n_peptides))
        return False

    # Meets all constraints
    logger.info("Configuration meets all constraints.")
    return True
