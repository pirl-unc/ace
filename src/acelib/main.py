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
The purpose of this python3 script is to implement the main API functions of ACE.
"""


import pandas as pd
import math
import random
import os
from golfy import init, is_valid, optimize
from ortools.sat.python import cp_model
from typing import Callable, List, Tuple
from .block_assignment import BlockAssignment
from .block_design import BlockDesign
from .constants import *
from .default_parameters import *
from .deconvolution import deconvolve_hit_peptides
from .logger import get_logger
from .sequence_features import AceNeuralEngine
from .types import *
from .utilities import *


logger = get_logger(__name__)


def run_ace_golfy(
        block_design: BlockDesign,
        random_seed: int = GENERATE_RANDOM_SEED,
        max_iters: int = GENERATE_GOLFY_MAX_ITERS,
        init_mode: str = GENERATE_GOLFY_INIT_MODE,
        verbose: bool = True
) -> BlockAssignment:
    """
    Generate a block assignment using golfy.

    Parameters
    ----------
    block_design        :   BlockDesign object.
    random_seed         :   Random seed.
    max_iters           :   Number of maximum iterations for golfy.
    init_mode           :   Init mode.
    verbose             :   If True, prints messages.

    Returns
    -------
    block_assignment    :   BlockAssignment object.
    """
    if verbose:
        logger.info('Started running golfy.')

    # Step 1. Set random seed
    random.seed(random_seed)

    # Step 2. Create a DataFrame of peptide IDs and index IDs
    df_peptides = convert_peptides_to_dataframe(peptides=block_design.peptides)
    df_peptides['peptide_index'] = list(range(0, len(df_peptides)))
    preferred_neighbors = [] # list of tuples (peptide ID index, peptide ID index)
    for peptide_id_1, peptide_id_2 in block_design.preferred_peptide_pairs:
        peptide_index_1 = df_peptides.loc[df_peptides['peptide_id'] == peptide_id_1,'peptide_index'].values[0]
        peptide_index_2 = df_peptides.loc[df_peptides['peptide_id'] == peptide_id_2,'peptide_index'].values[0]
        preferred_neighbors.append((peptide_index_1, peptide_index_2))

    # Step 3. Run golfy
    golfy_solution = init(
        num_peptides=len(df_peptides['peptide_id'].unique()),
        max_peptides_per_pool=block_design.num_peptides_per_pool,
        num_replicates=block_design.num_coverage,
        strategy=init_mode,
        preferred_neighbors=preferred_neighbors,
        verbose=verbose
    )
    optimize(golfy_solution, max_iters=max_iters, verbose=verbose)

    if verbose:
        logger.info('Finished running golfy.')

    # Step 4. Convert golfy assignments to a BlockAssignment object
    block_assignment = BlockAssignment.load_golfy_assignment(
        golfy_assignment=golfy_solution.assignments,
        df_peptides=df_peptides
    )

    return block_assignment


def run_ace_sat_solver(
        block_design: BlockDesign,
        max_peptides_per_pool: int = GENERATE_MAX_PEPTIDES_PER_POOL,
        num_processes: int = GENERATE_NUM_PROCESSES,
        shuffle_iters: int = GENERATE_SHUFFLE_ITERS,
        verbose: bool = True
) -> BlockAssignment:
    """
    Generate a block assignment using SAT solver.

    Parameters
    ----------
    block_design            :   BlockDesign object.
    max_peptides_per_pool   :   Maximum number of peptides per pool.
    num_processes           :   Number of processes (default: 4).
    shuffle_iters           :   Number of iterations to shuffle pool IDs (default: 100).
    verbose                 :   If True, prints messages.

    Returns
    -------
    block_assignment        :   BlockAssignment object.
    """
    if verbose:
        logger.info('Started running SAT solver.')

    # Step 1. Generate first coverage block assignment if there are preferred peptide pairs
    if len(block_design.preferred_peptide_pairs) > 0:
        block_assignment_1x_coverage = BlockAssignment.generate_single_coverage_block_assignment(
            peptides=block_design.peptides,
            preferred_peptide_pairs=block_design.preferred_peptide_pairs,
            num_peptides_per_pool=block_design.num_peptides_per_pool,
            coverage=1
        )
        block_design.disallowed_peptide_pairs = block_assignment_1x_coverage.pooled_peptide_pairs
        block_design.num_coverage -= 1
        start_pool_num = block_assignment_1x_coverage.num_pools + 1
        start_coverage_num = 2
    else:
        block_assignment_1x_coverage = BlockAssignment()
        start_pool_num = 1
        start_coverage_num = 1

    # Step 2. Divide block design into smaller computationally tractable designs.
    block_designs = BlockDesign.divide_block_design(
        block_design=block_design,
        max_peptides_per_block=block_design.max_peptides_per_block,
        max_peptides_per_pool=max_peptides_per_pool,
        verbose=verbose
    )
    if len(block_designs) > 1:
        if verbose:
            logger.info("Requested block design (%i peptides, %i peptides per pool, "
                        "%ix coverage) has been divided into the following smaller designs:" %
                        (block_design.num_peptides,
                         block_design.num_peptides_per_pool,
                         block_design.num_coverage))
            for block_designs_ in block_designs:
                logger.info('\tBlock assignments from the following block designs will be merged into one:')
                for block_design in block_designs_:
                    logger.info("\t\t%i peptides, %i peptide(s) per pool, %ix coverage, %i max peptides per block" %
                                (block_design.num_peptides,
                                 block_design.num_peptides_per_pool,
                                 block_design.num_coverage,
                                 block_design.max_peptides_per_block))

    # Step 3. Generate block assignments
    block_assignments = [block_assignment_1x_coverage]
    for block_designs_ in block_designs:
        start_pool_num_ = start_pool_num
        for block_design in block_designs_:
            if verbose:
                logger.info('Generating assignment for '
                            '%i peptides, '
                            '%i peptides per pool, '
                            '%ix coverage (%i maximum peptides per block)' %
                            (block_design.num_peptides,
                             block_design.num_peptides_per_pool,
                             block_design.num_coverage,
                             block_design.max_peptides_per_block))
            block_assignment = block_design.generate(
                random_seed=generate_random_seed(),
                num_processes=num_processes,
                verbose=verbose
            )
            block_assignment.assignments = BlockAssignment.update_ids(
                assignments=block_assignment.assignments,
                start_pool_num=start_pool_num_,
                start_coverage_num=start_coverage_num
            )
            start_pool_num_ += block_assignment.num_pools
            block_assignments.append(block_assignment)

    # Step 4. Merge assignments
    block_assignments = BlockAssignment.minimize_violations(
        block_assignments=block_assignments,
        shuffle_iters=shuffle_iters,
        verbose=verbose
    )
    block_assignment = BlockAssignment.merge(block_assignments=block_assignments)

    if verbose:
        logger.info('Finished running SAT solver.')
        logger.info('The returning block assignment has %i pools in total.' % block_assignment.num_pools)
        logger.info('The returning block assignment has %i peptides in total.' % len(block_assignment.peptide_ids))
    return block_assignment

