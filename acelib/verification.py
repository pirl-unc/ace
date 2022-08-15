#!/usr/bin/python3

"""
The purpose of this python3 script is to implement functions related to
verifying configuration and identifying hit peptide IDs.

Author: Jin Seok (Andy) Lee, Dhuvarakesh Karthikeyan

Last updated date: Aug 15, 2022
"""


import pandas as pd


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


# def verify_configuration_integrity(df_configuration):
#     """
#
#     """