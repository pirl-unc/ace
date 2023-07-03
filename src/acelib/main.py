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
from ortools.sat.python import cp_model
from typing import Callable, List, Tuple
from .constants import PlateTypes
from .elispot import ELISpot
from .logger import get_logger


logger = get_logger(__name__)


def run_ace_generate(
        df_peptides: pd.DataFrame,
        num_peptides_per_pool: int,
        num_coverage: int,
        num_processes: int,
        random_seed: int,
        disallowed_peptide_pairs: List[Tuple[str, str]] = [],
        enforced_peptide_pairs: List[Tuple[str, str]] = [],
        assign_well_ids: bool = True,
        plate_type: str = PlateTypes.PLATE_96_WELLS
) -> Tuple[int, pd.DataFrame]:
    """
    Generate an ELISpot configuration.

    Parameters
    ----------
    df_peptides                     :   pd.DataFrame with the following columns:
                                        'peptide_id'
                                        'peptide_sequence'
    num_peptides_per_pool           :   Number of peptides per pool.
    num_coverage                    :   Number of coverage (i.e. number of peptide replicates).
    num_processes                   :   Number of processes.
    random_seed                     :   Random seed.
    disallowed_peptide_pairs        :   List of tuples (peptide ID, peptide ID).
    enforced_peptide_pairs          :   List of tuples (peptide ID, peptide ID).
    assign_well_ids                 :   If true, assigns plate and well IDs to each pool ID.
    plate_type                      :   Plate type (allowed values: '96-well plate').

    Returns
    -------
    status                          :   Solver status; one of the following:
                                        'cp_model.OPTIMAL'
                                        'cp_model.FEASIBLE'
                                        'cp_model.INFEASIBLE'
                                        'cp_model.MODEL_INVALID'
                                        'cp_model.UNKNOWN'
    df_configuration                :   pd.DataFrame with the following columns:
                                        'pool_id',
                                        'peptide_id'
    """
    # Step 1. Create an ELISpot configuration
    elispot = ELISpot(
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage,
        num_processes=num_processes,
        peptide_ids=list(df_peptides['peptide_id'].unique())
    )

    # Step 2. Generate an ELISpot experiment configuration
    while True:
        status, df_configuration = elispot.generate_configuration(
            random_seed=random_seed,
            disallowed_peptide_pairs=disallowed_peptide_pairs,
            enforced_peptide_pairs=enforced_peptide_pairs
        )
        if status == cp_model.OPTIMAL:
            logger.info('An optimal configuration has been generated.')
            break
        else:
            logger.info('An optimal configuration could not be generated.')
            if len(disallowed_peptide_pairs) == 0 and len(enforced_peptide_pairs) == 0:
                break

        if len(disallowed_peptide_pairs) > 0:
            logger.info('Removing the last element in the list of disallowed peptide pairs as a constraint.')
            disallowed_peptide_pairs.pop()
        if len(enforced_peptide_pairs) > 0:
            logger.info('Removing the last element in the list of enforced peptide pairs as a constraint.')
            enforced_peptide_pairs.pop()

    # Step 3. Assign plate and well IDs
    if status == cp_model.OPTIMAL:
        if assign_well_ids:
            df_configuration = ELISpot.assign_well_ids(
                df_configuration=df_configuration,
                plate_type=plate_type
            )

    return status, df_configuration


def run_ace_identify(
        hit_pool_ids: List[str],
        df_configuration: pd.DataFrame,
) -> pd.DataFrame:
    """
    Identifies hit peptide IDs.

    Parameters
    ----------
    hit_pool_ids                    :   Hit pool IDs.
    df_configuration                :   DataFrame with the following columns:
                                        'pool_id'
                                        'peptide_id'

    Returns
    -------
    df_hits_max                     :   DataFrame with the following columns:
                                        'peptide_id'
                                        'pool_ids'
                                        'num_coverage'
                                        'deconvolution_result'
    """
    return ELISpot.identify_hit_peptides(
        hit_pool_ids=hit_pool_ids,
        df_configuration=df_configuration
    )


def run_ace_verify(
        df_configuration: pd.DataFrame,
        num_peptides_per_pool: int,
        num_coverage: int
) -> bool:
    """
    Verifies whether an ELIspot configuration meets all ACE constraints.

    Parameters
    ----------
    df_configuration        :   DataFrame with the following columns:
                                'pool_id'
                                'peptide_id'
    num_peptides_per_pool   :   Number of peptides per pool.
    num_coverage            :   Number of coverage (i.e. number of replicates per peptide).

    Returns
    -------
    is_valid                :   True if meets all ACE constraints.
                                False otherwise.
    """
    return ELISpot.verify_configuration(
        df_configuration=df_configuration,
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage
    )
