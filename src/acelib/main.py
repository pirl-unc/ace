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


import pandas as pd
import math
import random
import os
from golfy import init, is_valid, optimize
from ortools.sat.python import cp_model
from typing import Callable, List, Tuple
from .constants import *
from .default_parameters import *
from .elispot import ELISpot
from .logger import get_logger
from .sequence_features import AceNeuralEngine
from .utilities import convert_golfy_results, split_peptides


logger = get_logger(__name__)


def run_ace_golfy(
        df_peptides: pd.DataFrame,
        num_peptides_per_pool: int,
        num_coverage: int,
        random_seed:int,
        max_iters: int,
        init_mode: str,
        preferred_peptide_pairs: List[Tuple[str, str, float]] = []
) -> Tuple[bool, pd.DataFrame]:
    """
    Generate an ELISpot configuration using golfy.

    Parameters
    ----------
    df_peptides                 :   pd.DataFrame with the following columns:
                                    'peptide_id'
                                    'peptide_sequence'
    num_peptides_per_pool       :   Number of peptides per pool.
    num_coverage                :   Number of coverage (i.e. number of peptide replicates).
    random_seed                 :   Random seed.
    max_iters                   :   Number of maximum iterations for golfy.
    init_mode                   :   Init mode.
    preferred_peptide_pairs     :   List of tuples (peptide ID, peptide ID, score).

    Returns
    -------
    is_valid                    :   True if all peptides are put into a unique
                                    combination of pool IDs. False otherwise.
    df_configuration            :   pd.DataFrame with the following columns:
                                    'coverage_id'
                                    'pool_id'
                                    'peptide_id'
    """
    logger.info('Started running golfy.')

    # Step 1. Set random seed
    random.seed(random_seed)

    # Step 2. Create a DataFrame of peptide IDs and index IDs
    df_peptides['peptide_index'] = list(range(0, len(df_peptides)))
    preferred_neighbors = [] # list of tuples (peptide ID index, peptide ID index)
    for peptide_id_1, peptide_id_2, score in preferred_peptide_pairs:
        peptide_index_1 = df_peptides.loc[df_peptides['peptide_id'] == peptide_id_1,'peptide_index'].values[0]
        peptide_index_2 = df_peptides.loc[df_peptides['peptide_id'] == peptide_id_2,'peptide_index'].values[0]
        preferred_neighbors.append((peptide_index_1, peptide_index_2))

    # Step 3. Run golfy
    golfy_solution = init(
        num_peptides=len(df_peptides['peptide_id'].unique()),
        peptides_per_pool=num_peptides_per_pool,
        num_replicates=num_coverage,
        strategy=init_mode,
        preferred_neighbors=preferred_neighbors
    )
    optimize(golfy_solution, max_iters=max_iters)

    # Step 4. Convert golfy assignments to a DataFrame
    df_configuration = convert_golfy_results(
        golfy_assignment=golfy_solution.assignments,
        df_peptides=df_peptides
    )

    logger.info('Finished running golfy.')
    return is_valid(golfy_solution), df_configuration


def __run_sat_solver(
        df_peptides: pd.DataFrame,
        num_peptides_per_pool: int,
        num_coverage: int,
        num_processes: int,
        random_seed: int,
        is_first_coverage: bool,
        disallowed_peptide_pairs: List[Tuple[str, str]] = []
) -> Tuple[int, pd.DataFrame]:
    """
    Generates an ELISpot configuration by running a SAT (social golfer problem) solver.

    Parameters
    ----------
    df_peptides                     :   pd.DataFrame with the following columns:
                                        'peptide_id'
                                        'peptide_sequence'
    num_peptides_per_pool           :   Number of peptides per pool.
    num_coverage                    :   Number of coverage (i.e. number of peptide replicates).
    num_processes                   :   Number of processes.
    random_seed                     :   Random seed.
    is_first_coverage               :   True if first coverage. False otherwise.
    disallowed_peptide_pairs        :   List of tuples (peptide ID, peptide ID).

    Returns
    -------
    status                          :   Solver status; one of the following:
                                        'cp_model.OPTIMAL'
                                        'cp_model.FEASIBLE'
                                        'cp_model.INFEASIBLE'
                                        'cp_model.MODEL_INVALID'
                                        'cp_model.UNKNOWN'
    df_configuration                :   pd.DataFrame with the following columns:
                                        'coverage_id'
                                        'pool_id',
                                        'peptide_id'
                                        'peptide_sequence'
    """
    # Step 1. Create an ELISpot configuration
    elispot = ELISpot(
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage,
        num_processes=num_processes,
        df_peptides=df_peptides
    )

    # Step 2. Generate an ELISpot experiment configuration
    status, df_configuration = elispot.generate_configuration(
        random_seed=random_seed,
        is_first_coverage=is_first_coverage,
        disallowed_peptide_pairs=disallowed_peptide_pairs
    )
    return status, df_configuration


def run_ace_sat_solver(
        df_peptides: pd.DataFrame,
        num_peptides_per_pool: int,
        num_coverage: int,
        num_peptides_per_batch: int,
        random_seed: int,
        num_processes: int,
        preferred_peptide_pairs: List[Tuple[str,str,float]] = []
) -> pd.DataFrame:
    """
    Generate an ELISpot configuration using SAT solver.

    Parameters
    ----------
    df_peptides                 :   pd.DataFrame with the following columns:
                                    'peptide_id'
                                    'peptide_sequence'
    num_peptides_per_pool       :   Number of peptides per pool.
    num_coverage                :   Number of coverage (i.e. number of peptide replicates).
    num_peptides_per_batch      :   Number of peptides per batch.
    random_seed                 :   Random seed.
    num_processes               :   Number of processes.
    preferred_peptide_pairs     :   List of tuples (peptide ID, peptide ID, score).

    Returns
    -------
    df_configuration            :   pd.DataFrame with the following columns:
                                    'coverage_id'
                                    'pool_id'
                                    'peptide_id'
    """
    logger.info('Started running SAT solver.')

    # Generate first coverage configuration if there are preferred peptide pairs
    if len(preferred_peptide_pairs) > 0:
        df_configuration_first_coverage = ELISpot.generate_first_coverage_configuration(
            df_peptides=df_peptides,
            preferred_peptide_pairs=preferred_peptide_pairs,
            num_peptides_per_pool=num_peptides_per_pool
        )
        disallowed_peptide_pairs = ELISpot.fetch_pooled_peptide_pairs(
            df_configuration=df_configuration_first_coverage
        )
        is_first_coverage = False
        num_coverage = num_coverage - 1
    else:
        df_configuration_first_coverage = pd.DataFrame()
        disallowed_peptide_pairs = []
        is_first_coverage = True
        num_coverage = num_coverage

    if len(list(df_peptides['peptide_id'].unique())) > num_peptides_per_batch:
        # Split peptides into batches
        list_df = split_peptides(
            df_peptides=df_peptides,
            num_peptides_per_batch=num_peptides_per_batch
        )
        df_configuration = pd.DataFrame()
        for df_peptides_ in list_df:
            status, df_configuration_ = __run_sat_solver(
                df_peptides=df_peptides_,
                num_peptides_per_pool=num_peptides_per_pool,
                num_coverage=num_coverage,
                num_processes=num_processes,
                random_seed=random_seed,
                is_first_coverage=is_first_coverage,
                disallowed_peptide_pairs=disallowed_peptide_pairs
            )
            if status != cp_model.OPTIMAL:
                logger.error('Exiting program. Please review your configuration parameters before running SAT-solver again.')
                exit(1)
            df_configuration = pd.concat([df_configuration, df_configuration_])
    else:
        status, df_configuration = __run_sat_solver(
            df_peptides=df_peptides,
            num_peptides_per_pool=num_peptides_per_pool,
            num_coverage=num_coverage,
            num_processes=num_processes,
            random_seed=random_seed,
            is_first_coverage=is_first_coverage,
            disallowed_peptide_pairs=disallowed_peptide_pairs
        )
        if status != cp_model.OPTIMAL:
            logger.error('Exiting program. Please review your configuration parameters before running SAT-solver again.')
            exit(1)

    df_configuration = pd.concat([df_configuration_first_coverage, df_configuration])

    logger.info('Finished running SAT solver.')
    return df_configuration


def run_ace_deconvolve(
        hit_pool_ids: List[str],
        df_configuration: pd.DataFrame,
        min_coverage: int
) -> pd.DataFrame:
    """
    Deconvolves hit peptide IDs.

    Parameters
    ----------
    hit_pool_ids                    :   Hit pool IDs.
    df_configuration                :   DataFrame with the following columns:
                                        'coverage_id'
                                        'pool_id'
                                        'peptide_id'
                                        'peptide_sequence'
    min_coverage                    :   Minimum coverage.

    Returns
    -------
    df_hits                         :   DataFrame with the following columns:
                                        'peptide_id'
                                        'peptide_sequence'
                                        'pool_ids'
                                        'num_coverage'
                                        'deconvolution_result'
    """
    return ELISpot.deconvolve_hit_peptides(
        hit_pool_ids=hit_pool_ids,
        df_configuration=df_configuration,
        min_coverage=min_coverage
    )


def run_ace_verify(
        df_configuration: pd.DataFrame,
        num_peptides_per_pool: int,
        num_coverage: int
) -> bool:
    """
    Verifies whether an ELISpot configuration meets all criteria for an
    optimal configuration.

    Parameters
    ----------
    df_configuration        :   DataFrame with the following columns:
                                'pool_id'
                                'peptide_id'
    num_peptides_per_pool   :   Number of peptides per pool.
    num_coverage            :   Number of coverage (i.e. number of replicates per peptide).

    Returns
    -------
    is_optimal              :   True if meets all ACE constraints.
                                False otherwise.
    """
    return ELISpot.verify_configuration(
        df_configuration=df_configuration,
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage
    )

