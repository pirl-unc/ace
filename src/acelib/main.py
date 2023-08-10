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
import torch
from golfy import init, is_valid, optimize, deconvolve
from ortools.sat.python import cp_model
from transformers import AutoTokenizer, AutoModelForMaskedLM
from typing import Callable, List, Tuple
from .block_assignment import BlockAssignment
from .block_design import BlockDesign
from .constants import *
from .default_parameters import *
from .deconvolution import convert_to_golfy_spotcounts, empirical_deconvolve, DeconvolutionResult
from .logger import get_logger
from .sequence_features import AceNeuralEngine
from .types import *
from .utilities import *


logger = get_logger(__name__)


def run_ace_golfy(
        block_design: BlockDesign,
        random_seed: int = GENERATE_GOLFY_RANDOM_SEED,
        max_iters: int = GENERATE_GOLFY_MAX_ITERS,
        strategy: str = GENERATE_GOLFY_STRATEGY,
        allow_extra_pools: bool = GENERATE_GOLFY_ALLOW_EXTRA_POOLS,
        verbose: bool = True
) -> BlockAssignment:
    """
    Generate a block assignment using golfy.

    Parameters
    ----------
    block_design        :   BlockDesign object.
    random_seed         :   Random seed.
    max_iters           :   Number of maximum iterations for golfy (default: 2000).
    strategy            :   Initialization strategy (default: 'greedy').
    allow_extra_pools   :   Allow extra pools (default: False).
    verbose             :   Print logs (default: True).

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
        strategy=strategy,
        preferred_neighbors=preferred_neighbors,
        allow_extra_pools=allow_extra_pools,
        verbose=verbose
    )
    optimize(
        golfy_solution,
        max_iters=max_iters,
        allow_extra_pools=allow_extra_pools,
        verbose=verbose
    )

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
        max_peptides_per_pool: int = GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_POOL,
        num_processes: int = GENERATE_CPSAT_SOLVER_NUM_PROCESSES,
        shuffle_iters: int = GENERATE_CPSAT_SOLVER_SHUFFLE_ITERS,
        verbose: bool = True
) -> BlockAssignment:
    """
    Generate a block assignment using SAT solver.

    Parameters
    ----------
    block_design            :   BlockDesign object.
    max_peptides_per_pool   :   Maximum number of peptides per pool (default: 10).
    num_processes           :   Number of processes (default: 2).
    shuffle_iters           :   Number of iterations to shuffle pool IDs (default: 1000).
    verbose                 :   Print log (default: True).

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
    if verbose:
        logger.info('Started minimizing violations.')
    block_assignments = BlockAssignment.minimize_violations(
        block_assignments=block_assignments,
        shuffle_iters=shuffle_iters,
        verbose=verbose
    )
    if verbose:
        logger.info('Finished minimizing violations.')
    block_assignment = BlockAssignment.merge(block_assignments=block_assignments)

    if verbose:
        logger.info('Finished running SAT solver.')
        logger.info('The returning block assignment has the following:')
        logger.info('\t%i pools in total.' % block_assignment.num_pools)
        logger.info('\t%i peptides in total.' % len(block_assignment.peptide_ids))
    return block_assignment


def run_ace_generate(
        peptides: Peptides,
        num_peptides_per_pool: int,
        num_coverage: int,
        trained_model_file: str,
        cluster_peptides: bool,
        mode: GenerateModes.ALL = GenerateModes.GOLFY,
        sequence_similarity_function: SequenceSimilarityFunctions.ALL = SequenceSimilarityFunctions.EUCLIDEAN,
        sequence_similarity_threshold: float = 0.8,
        golfy_random_seed: int = generate_random_seed(),
        golfy_strategy: GolfyStrategies.ALL = GENERATE_GOLFY_STRATEGY,
        golfy_max_iters: int = GENERATE_GOLFY_MAX_ITERS,
        golfy_allow_extra_pools: bool = GENERATE_GOLFY_ALLOW_EXTRA_POOLS,
        cpsat_solver_num_processes: int = GENERATE_CPSAT_SOLVER_NUM_PROCESSES,
        cpsat_solver_shuffle_iters: int = GENERATE_CPSAT_SOLVER_SHUFFLE_ITERS,
        cpsat_solver_max_peptides_per_block: int = GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_BLOCK,
        cpsat_solver_max_peptides_per_pool: int = GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_POOL,
        assign_well_ids: bool = GENERATE_ASSIGN_WELL_IDS,
        num_plate_wells: PlateWells.ALL = PlateWells.WELLS_96, 
        verbose: bool = True
) -> Tuple[BlockAssignment, BlockDesign]:
    """
    Runs ACE 'generate' command.

    Parameters
    ----------
    peptides                            :   List of tuples (peptide ID, peptide sequence).
    num_peptides_per_pool               :   Number of peptides per pool.
    num_coverage                        :   Coverage.
    trained_model_file                  :   Trained model file.
    cluster_peptides                    :   Cluster peptides.
    mode                                :   'golfy' or 'cpsat_solver' (default: 'golfy').
    sequence_similarity_function        :   'euclidean', 'cosine' or 'levenshtein' (default: 'euclidean').
    sequence_similarity_threshold       :   Sequence similarity threshold. Recommended values: 
                                            0.8 for 'euclidean'
                                            0.9 for 'cosine'
                                            3 for 'levenshtein'
    golfy_random_seed                   :   Random seed for golfy (default: randomly generated).
    golfy_strategy                      :   'greedy', 'random', 'valid', 'singleton' or 'repeat' (default: 'greedy').
    golfy_max_iters                     :   Maximum number of iterations for golfy (default: 2000).
    golfy_allow_extra_pools             :   Allow extra pools for golfy (default: False).
    cpsat_solver_num_processes          :   Number of processes for CP-SAT solver (default: 2).
    cpsat_solver_shuffle_iters          :   Number of iterations to shuffle for CP-SAT solver (default: 1000).
    cpsat_solver_max_peptides_per_block :   Maximum number of peptides per block for CP-SAT solver (default: 10).
    cpsat_solver_max_peptides_per_pool  :   Maximum number of peptides per pool for CP-SAT solver (default: 100).
    assign_well_ids                     :   Assigns well IDs (default: True).
    num_plate_wells                     :   Number of wells on plate (default: 96).
    verbose                             :   Print logs (default: True).
                                        
    Returns
    -------
    block_assignment                    :   BlockAssignment object.
    block_design                        :   BlockDesign object.
    """
    # Step 2. Identify pairs of similar peptides
    if cluster_peptides:
        if sequence_similarity_function == SequenceSimilarityFunctions.LEVENSHTEIN:
            preferred_peptide_pairs = AceNeuralEngine.find_levenshtein_paired_peptides(
                peptide_ids=[p[0] for p in peptides],
                peptide_sequences=[p[1] for p in peptides],
                threshold=sequence_similarity_threshold
            )
        else:
            ESM2_TOKENIZER = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
            ESM2_MODEL = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D", return_dict=True, output_hidden_states=True)
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            ace_eng = AceNeuralEngine(ESM2_MODEL, ESM2_TOKENIZER, device)
            ace_eng.load_weights(trained_model_file)
            preferred_peptide_pairs = ace_eng.find_paired_peptides(
                peptide_ids=[p[0] for p in peptides],
                peptide_sequences=[p[1] for p in peptides],
                sim_fxn=sequence_similarity_function,
                threshold=sequence_similarity_threshold
            )
        if verbose:
            logger.info('%i peptide cluster(s) identified by the ACE sequence similarity neural engine:' % len(preferred_peptide_pairs))
            logger.info('peptide ID, peptide ID: similarity score')
            for peptide_id_1, peptide_id_2, score in preferred_peptide_pairs:
                logger.info('%s, %s: %f' % (peptide_id_1, peptide_id_2, score))
    else:
        preferred_peptide_pairs = []

    preferred_peptide_pairs = [(p1, p2) for p1, p2, score in preferred_peptide_pairs]

    # Step 3. Generate a block design
    block_design = BlockDesign(
        peptides=peptides,
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage,
        max_peptides_per_block=cpsat_solver_max_peptides_per_block,
        disallowed_peptide_pairs=[],
        allow_extra_pools=golfy_allow_extra_pools,
        preferred_peptide_pairs=preferred_peptide_pairs
    )

    # Step 4. Generate a block assignment
    if mode == GenerateModes.GOLFY:
        block_assignment = run_ace_golfy(
            block_design=block_design,
            random_seed=golfy_random_seed,
            max_iters=golfy_max_iters,
            strategy=golfy_strategy,
            allow_extra_pools=golfy_allow_extra_pools,
            verbose=verbose
        )
    elif mode == GenerateModes.CPSAT_SOLVER:
        block_assignment = run_ace_sat_solver(
            block_design=block_design,
            max_peptides_per_pool=cpsat_solver_max_peptides_per_pool,
            num_processes=cpsat_solver_num_processes,
            shuffle_iters=cpsat_solver_shuffle_iters,
            verbose=verbose
        )
    else:
        logger.error('Unknown mode: %s' % mode)
        exit(1)

    # Step 5. Check if the block assignment is optimal
    # re-assign value because it could have been decremented
    # to account for preferred peptide neighbors (1x coverage).
    block_design.num_coverage = num_coverage
    block_assignment.is_optimal(
        num_coverage=block_design.num_coverage,
        num_peptides_per_pool=block_design.num_peptides_per_pool,
        verbose=verbose
    )

    # Step 6. Assign plate and well IDs
    if assign_well_ids:
        block_assignment.assign_well_ids(num_plate_wells=num_plate_wells)
    
    return block_assignment, block_design


def run_ace_deconvolve(
        df_readout: pd.DataFrame,
        block_assignment: BlockAssignment,
        mode,
        statistical_min_peptide_activity: float,
        empirical_min_coverage: int,
        empirical_min_spot_count: int,
        verbose: bool = True
) -> DeconvolutionResult:
    """
    Deconvolves pool spot counts.

    Parameters
    ----------
    df_readout      :   pd.DataFrame with the following columns:
                        'pool_id'
                        'spot_count'
    df_assignment   :   pd.DataFrame with the following columns:
                        'coverage_id'
                        'pool_id'
                        'peptide_id'

    Returns
    -------

    """
    if mode == DeconvolveModes.EMPIRICAL:
        df_readout = df_readout[df_readout['spot_count'] >= empirical_min_spot_count]
        hit_pool_ids = list(df_readout['pool_id'].unique())
        return empirical_deconvolve(
            hit_pool_ids=hit_pool_ids,
            df_assignment=block_assignment.to_dataframe(),
            min_coverage=empirical_min_coverage
        )
    elif mode == DeconvolveModes.EM or mode == DeconvolveModes.LASSO:
        golfy_design, peptide_indices = block_assignment.to_golfy_design()
        spot_counts = {}
        for _, row in df_readout.iterrows():
            spot_counts[int(row['pool_id'])] = int(row['spot_count'])
        golfy_spot_counts = convert_to_golfy_spotcounts(
            spot_counts=spot_counts,
            block_assignment=block_assignment
        )
        golfy_deconvolve_result = deconvolve(
            s=golfy_design,
            spot_counts=golfy_spot_counts,
            method=mode,
            min_peptide_activity=statistical_min_peptide_activity,
            verbose=verbose
        )
        deconvolution_result = DeconvolutionResult()
        peptide_idx = 0
        for peptide_activity in golfy_deconvolve_result.activity_per_peptide:
            peptide_id = peptide_indices[peptide_idx]
            peptide_sequence = block_assignment.get_peptide_sequence(peptide_id=peptide_id)
            if peptide_idx in golfy_deconvolve_result.high_confidence_hits:
                label = DeconvolutionLabels.CONFIDENT_HIT
            else:
                label = DeconvolutionLabels.NOT_A_HIT
            deconvolution_result.add_peptide(
                peptide_id=peptide_id,
                peptide_sequence=peptide_sequence,
                peptide_activity_level=peptide_activity,
                label=label,
                pool_ids=[]
            )
            peptide_idx += 1
        return deconvolution_result
    else:
        logger.error('Unknown mode: %s' % mode)
        exit(1)

