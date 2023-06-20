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
from .elispot import ELIspot
from .logger import get_logger


logger = get_logger(__name__)


def run_ace_generate(
        df_peptides: pd.DataFrame,
        num_peptides_per_pool: int,
        num_coverage: int,
        num_processes: int,
        dissimilarity_inference_func: Callable[[pd.DataFrame], List[Tuple[str,str]]] = None,
        assign_well_ids: bool = True,
        plate_type: str = PlateTypes.PLATE_96_WELLS
) -> Tuple[int, pd.DataFrame]:
    """
    Generate an ELIspot configuration.

    Parameters
    ----------
    df_peptides                     :   pd.DataFrame with the following columns:
                                        'peptide_id'
                                        'peptide_sequence'
    num_peptides_per_pool           :   Number of peptides per pool.
    num_coverage                    :   Number of coverage (i.e. number of peptide replicates).
    num_processes                   :   Number of processes.
    dissimilarity_inference_func    :   Custom function to predicting dissimilarity between two peptides.
                                        The function takes as input a pandas DataFrame with the following columns:
                                        'peptide_id', 'peptide_'sequence'.
                                        The function outputs a list of tuples (peptide ID 1, peptide ID 2)
                                        in the descending order of similarity.
                                        If CP-SAT solver cannot find an optimal solution, the experiment
                                        generation will iteratively forgo the last element in the sorted list as a constraint.
    assign_well_ids                 :   If true, assigns plate and well IDs to each pool ID.
    plate_type                      :   Plate type (allowed values: '96-well plate').

    Returns
    -------
    status
    df_configuration        :   pd.DataFrame with the following columns:
                                'pool_id',
                                'peptide_id'
    """
    # Step 1. Create an ELIspot configuration
    elispot = ELIspot(
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage,
        num_processes=num_processes,
        peptide_ids=list(df_peptides['peptide_id'].unique())
    )

    # Step 2. Identify disallowed peptide pairs
    if dissimilarity_inference_func is not None:
        disallowed_peptide_pairs = dissimilarity_inference_func(df_peptides)
    else:
        disallowed_peptide_pairs = []

    # Step 3. Generate an ELIspot experiment configuration
    while True:
        status, df_configuration = elispot.generate_configuration(
            disallowed_peptide_pairs=disallowed_peptide_pairs
        )
        if status == cp_model.OPTIMAL:
            logger.info('An optimal configuration has been generated.')
            break
        if len(disallowed_peptide_pairs) > 0:
            logger.info('An optimal configuration could not be generated. Removing the last element in disallowed peptide pairs (i.e. least similar) as a constraint.')
            disallowed_peptide_pairs.pop()
        else:
            logger.info('An optimal configuration could not be generated.')
            break

    # Step 4. Assign plate and well IDs
    if status == cp_model.OPTIMAL:
        if assign_well_ids:
            df_configuration = ELIspot.assign_well_ids(
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
    return ELIspot.identify_hit_peptides(
        hit_pool_ids=hit_pool_ids,
        df_configuration=df_configuration
    )


def run_ace_generate_with_precomputed_configuration(
        df_peptides: pd.DataFrame,
        df_template_configuration: pd.DataFrame,
        num_peptides_per_pool: int,
        num_coverage: int,
        dissimilarity_inference_func: Callable[[pd.DataFrame], List[Tuple[str,str]]] = None,
        assign_well_ids: bool = True,
        plate_type: str = PlateTypes.PLATE_96_WELLS
) -> Tuple[int, pd.DataFrame]:
    """
    Generate an ELIspot configuration by recycling a pre-computed configuration.

    Parameters
    ----------
    df_peptides                     :   DataFrame with the following columns:
                                        'peptide_id'
                                        'peptide_sequence'
    df_template_configuration       :   DataFrame with the following columns:
                                        'coverage_id'
                                        'pool_id'
                                        'peptide_id'
    num_peptides_per_pool           :   Number of peptides per pool.
    num_coverage                    :   Number of coverage (i.e. number of peptide replicates).
    dissimilarity_inference_func    :   Custom function to predicting dissimilarity between two peptides.
                                        The function takes as input a pandas DataFrame with the following columns:
                                        'peptide_id', 'peptide_'sequence'.
                                        The function outputs a list of tuples (peptide ID 1, peptide ID 2)
                                        in the descending order of similarity.
                                        If CP-SAT solver cannot find an optimal solution, the experiment
                                        generation will iteratively forgo the last element in the sorted list as a constraint.
    assign_well_ids                 :   If true, assigns plate and well IDs to each pool ID.
    plate_type                      :   Plate type (allowed values: '96-well plate').

    Returns
    -------
    status
    df_configuration        :   pd.DataFrame with the following columns:
                                'pool_id',
                                'peptide_id'
    """
    # Step 1. Create an ELIspot configuration
    elispot = ELIspot(
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage,
        num_processes=1,
        peptide_ids=list(df_peptides['peptide_id'].unique())
    )

    # Step 2. Identify disallowed peptide pairs
    if dissimilarity_inference_func is not None:
        disallowed_peptide_pairs = dissimilarity_inference_func(df_peptides)
    else:
        disallowed_peptide_pairs = []

    # Step 3. Generate an ELIspot experiment configuration
    df_configuration = elispot.recycle_configuration(
        df_template_configuration=df_template_configuration,
        disallowed_peptide_pairs=disallowed_peptide_pairs
    )

    # Step 4. Assign plate and well IDs
    if assign_well_ids:
        df_configuration = ELIspot.assign_well_ids(
            df_configuration=df_configuration,
            plate_type=plate_type
        )
    return df_configuration


def run_ace_verify(
        df_configuration: pd.DataFrame,
        num_coverage: int
) -> bool:
    """
    Verifies whether an ELIspot configuration meets all ACE constraints.

    Parameters
    ----------
    df_configuration    :   DataFrame with the following columns:
                            'pool_id'
                            'peptide_id'
    num_coverage        :   Number of coverage (i.e. number of replicates per peptide).

    Returns
    -------
    is_valid            :   True if meets all ACE constraints.
                            False otherwise.
    """
    return ELIspot.verify_configuration(
        df_configuration=df_configuration,
        num_coverage=num_coverage
    )
