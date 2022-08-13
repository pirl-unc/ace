#!/usr/bin/python3

"""
The purpose of this python3 script is to implement functions related to
solving the constraint satisfaction problem on pool membership assignment.

Author: Jin Seok (Andy) Lee, Dhuvarakesh Karthikeyan

Last updated date: Aug 12, 2022
"""


import pandas as pd
from ortools.sat.python import cp_model
from itertools import product, groupby, combinations
from acelib.logger import get_logger


logger = get_logger(__name__)


def solve_assay_configuration(n_peptides: int,
                              n_peptides_per_pool: int,
                              n_coverage: int,
                              peptide_ids: list = []):
    """
    This functions solves the assay configuration as a constraint problem.

    Args
    ----
    n_peptides              :   Number of total peptides.
    n_peptides_per_pool     :   Number of peptides per pool.
                                Make sure this integer is a factor of n_peptides.
    n_coverage              :   Coverage.
    peptide_ids             :   List of peptide IDs (optional). If this is
                                specified, 'n_peptides' does not have to be specified.

    Returns
    -------
    DataFrame with the following columns:
    'coverage', 'pool_id', 'peptide_id'
    """
    # Step 1. Check if peptide IDs have been specified.
    if len(peptide_ids) > 0:
        n_peptides = len(peptide_ids)

    # Step 2. Make sure 'n_peptides_per_pool' is a factor of 'n_peptides'
    if n_peptides % n_peptides_per_pool != 0:
        logger.error("The parameter 'n_peptides_per_pool' must be a factor of 'n_peptides'.")
        return

    # Step 3. Derive the number of pools per coverage
    n_pools = int(n_peptides / n_peptides_per_pool) # per coverage
    logger.info("Number of peptides: %i." % n_peptides)
    logger.info("Number of coverage: %i." % n_coverage)
    logger.info("Number of peptides per pool: %i." % n_peptides_per_pool)
    logger.info("Number of pools per coverage: %i." % n_pools)

    # Step 4. Generate IDs for peptides, pools, and coverage
    if len(peptide_ids) == 0:
        peptide_ids = ['peptide_' + str(i) for i in range(1, n_peptides + 1)]
    pool_ids = ['pool_' + str(i) for i in range(1, n_pools + 1)]
    coverage_ids = ['coverage_' + str(i) for i in range(1, n_coverage + 1)]

    # Step 5. Construct a constraint programming model
    model = cp_model.CpModel()
    data_dict = {
        'combination_id': [],
        'coverage_id': [],
        'pool_id': [],
        'peptide_id': [],
        'bool_variable': []
    }
    var_dict = {}
    for curr_coverage_id in coverage_ids:
        for curr_pool_id in pool_ids:
            for curr_peptide_id in peptide_ids:
                curr_combination_id = "%s/%s/%s" % (curr_coverage_id, curr_pool_id, curr_peptide_id)
                curr_bool_var = model.NewBoolVar(curr_combination_id)
                data_dict['combination_id'].append(curr_combination_id)
                data_dict['coverage_id'].append(curr_coverage_id)
                data_dict['pool_id'].append(curr_pool_id)
                data_dict['peptide_id'].append(curr_peptide_id)
                data_dict['bool_variable'].append(curr_bool_var)
                var_dict[(curr_coverage_id, curr_pool_id, curr_peptide_id)] = curr_bool_var
    df = pd.DataFrame(data_dict)

    # Constraint 1: each peptide is in only one pool in each coverage
    for name, group in df.groupby(['coverage_id', 'peptide_id']):
        model.Add(sum(row['bool_variable'] for index, row in group.iterrows()) == 1)

    # Constraint 2: each pool consists of the fixed number of peptides
    for name, group in df.groupby(['coverage_id', 'pool_id']):
        model.Add(sum(row['bool_variable'] for index, row in group.iterrows()) == n_peptides_per_pool)

    # Constraint 3: two peptides cannot be in the same pool more than once across all coverages
    for peptide_id_1, peptide_id_2 in combinations(peptide_ids, r=2):
        peptide_pair_bool_variables = []
        for curr_coverage_id in coverage_ids:
            for curr_pool_id in pool_ids:
                pair_bool_variable = model.NewBoolVar("%s/%s/%s/%s" % (curr_coverage_id, curr_pool_id, peptide_id_1, peptide_id_2))
                peptide_1_bool_variable = var_dict[(curr_coverage_id, curr_pool_id, peptide_id_1)]
                peptide_2_bool_variable = var_dict[(curr_coverage_id, curr_pool_id, peptide_id_2)]

                # pair_bool_variable has to be 1 if peptide 1 and peptide 2 are paired together
                model.Add((peptide_1_bool_variable + peptide_2_bool_variable - pair_bool_variable) <= 1)
                peptide_pair_bool_variables.append(pair_bool_variable)
        # All pairs can appear once
        model.Add(sum(peptide_pair_bool_variables) <= 1)

    # Step 6. Solve
    logger.info("CP solver started")
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 8
    solver.Solve(model)
    # print(solver.ResponseStats())
    logger.info("CP solver finished")

    # Step 7. Parse solution
    solutions_data = {
        'coverage_id': [],
        'pool_id': [],
        'peptide_id': []
    }
    for curr_bool_variable in df['bool_variable'].values.tolist():
        if solver.Value(curr_bool_variable) == 1:
            curr_bool_variable_elements = str(curr_bool_variable).split("/")
            curr_coverage_id = curr_bool_variable_elements[0]
            curr_pool_id = curr_bool_variable_elements[1]
            curr_peptide_id = curr_bool_variable_elements[2]
            solutions_data['coverage_id'].append(curr_coverage_id)
            solutions_data['pool_id'].append(curr_pool_id)
            solutions_data['peptide_id'].append(curr_peptide_id)
    df_solutions = pd.DataFrame(solutions_data)
    return df_solutions
