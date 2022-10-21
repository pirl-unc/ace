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
The purpose of this python3 script is to implement functions related to
solving the constraint satisfaction problem on pool membership assignment.
"""


from __future__ import print_function, division, absolute_import


import pandas as pd
from itertools import combinations
from ortools.sat.python import cp_model
from .logging import get_logger
from .default_parameters import *


logger = get_logger(__name__)


def generate_assay_configuration(n_peptides: int,
                                 n_peptides_per_pool: int,
                                 n_coverage: int,
                                 peptide_ids: list = [],
                                 disallowed_peptide_pairs: list = [],
                                 num_processes=NUM_PROCESSES):
    """
    Generates an assay configuration.

    Parameters
    ----------
    n_peptides               :  Number of total peptides.
    n_peptides_per_pool      :  Number of peptides per pool.
                                Make sure this integer is a factor of n_peptides.
    n_coverage               :  Coverage.
    peptide_ids              :  List of peptide IDs (optional).
                                If this is specified, 'n_peptides'
                                does not have to be specified.
    disallowed_peptide_pairs :  Disallowed peptide pairs (optional).
                                Each value is a tuple
                                (e.g. [('peptide_1', 'peptide_2')]).
    num_processes            :  Number of processes.

    Returns
    -------
    df_configuration         :  DataFrame with the following columns:
                                'coverage_id'
                                'pool_id'
                                'peptide_id'
    """
    # Step 1. Check if peptide IDs have been specified.
    if len(peptide_ids) > 0:
        n_peptides = len(peptide_ids)

    # Step 2. Make sure 'n_peptides_per_pool' is a factor of 'n_peptides'
    if n_peptides % n_peptides_per_pool != 0:
        logger.error("The parameter 'n_peptides_per_pool' must be a factor of 'n_peptides'.")
        return

    # Step 3. Derive the number of pools per coverage
    n_pools = int(n_peptides / n_peptides_per_pool)
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
        'coverage_id': [],
        'pool_id': [],
        'peptide_id': [],
        'bool_variable': []
    }
    # Initialize the constraint programming model dictionary
    var_dict = {}
    for curr_coverage_id in coverage_ids:
        for curr_pool_id in pool_ids:
            for curr_peptide_id in peptide_ids:
                curr_bool_var = model.NewBoolVar("%s/%s/%s" % (curr_coverage_id, curr_pool_id, curr_peptide_id))
                data_dict['coverage_id'].append(curr_coverage_id)
                data_dict['pool_id'].append(curr_pool_id)
                data_dict['peptide_id'].append(curr_peptide_id)
                data_dict['bool_variable'].append(curr_bool_var)
                var_dict[(curr_coverage_id, curr_pool_id, curr_peptide_id)] = curr_bool_var
    df = pd.DataFrame(data_dict)

    # Constraint 1. Each peptide appears exactly once in each coverage
    for curr_coverage_id in df['coverage_id'].unique():
        for curr_peptide_id in df['peptide_id'].unique():
            df_matched = df.loc[(df['coverage_id'] == curr_coverage_id) &
                                (df['peptide_id'] == curr_peptide_id),:]
            model.Add(sum(df_matched['bool_variable'].values.tolist()) == 1)

    # Constraint 2. Each pool in each coverage has exactly the same number of peptides
    for curr_coverage_id in df['coverage_id'].unique():
        for curr_pool_id in df['pool_id'].unique():
            df_matched = df.loc[(df['coverage_id'] == curr_coverage_id) &
                                (df['pool_id'] == curr_pool_id),:]
            model.Add(sum(df_matched['bool_variable'].values.tolist()) == n_peptides_per_pool)

    # Constraint 3. No two peptides are in the same pool more than once
    for peptide_id_1, peptide_id_2 in combinations(peptide_ids, r=2):
        peptide_pair_bool_variables = []
        for curr_coverage_id in coverage_ids:
            for curr_pool_id in pool_ids:
                pair_bool_variable = model.NewBoolVar("%s/%s/%s" % (curr_coverage_id, peptide_id_1, peptide_id_2))
                peptide_1_bool_variable = var_dict[(curr_coverage_id, curr_pool_id, peptide_id_1)]
                peptide_2_bool_variable = var_dict[(curr_coverage_id, curr_pool_id, peptide_id_2)]

                # pair_bool_variable has to be 1 if peptide 1 and peptide 2 are paired together
                model.Add((peptide_1_bool_variable + peptide_2_bool_variable - pair_bool_variable) <= 1)
                peptide_pair_bool_variables.append(pair_bool_variable)
        # All pairs can appear together in the same pool at most once
        model.Add(sum(peptide_pair_bool_variables) <= 1)

    # Constraint 4. Apply constraints for disallowed peptide pairs



    # Step 6. Solve
    logger.info("CP solver started")
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = num_processes
    solver.enumerate_all_solutions = False
    status = solver.Solve(model)
    logger.info("CP solver finished")

    if status == cp_model.OPTIMAL:
        logger.info("Solution is optimal.")
    else:
        logger.info("Solution is not optimal.")

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

            # Fix pool ID
            if curr_coverage_id != "coverage_1":
                curr_coverage_id_int = int(curr_coverage_id.split('_')[1])
                curr_pool_id = 'pool_' + str(n_pools * (curr_coverage_id_int - 1) + int(curr_pool_id.split('_')[1]))

            solutions_data['coverage_id'].append(curr_coverage_id)
            solutions_data['pool_id'].append(curr_pool_id)
            solutions_data['peptide_id'].append(curr_peptide_id)
    df_configuration = pd.DataFrame(solutions_data)
    return df_configuration

