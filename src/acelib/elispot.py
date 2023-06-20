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
The purpose of this python3 script is to implement the ELIspot dataclass.
"""


import itertools
import math
import pandas as pd
import random
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations, product
from ortools.sat.python import cp_model
from typing import List, Tuple
from .constants import DeconvolutionResults, PlateTypes
from .logger import get_logger


logger = get_logger(__name__)


@dataclass(frozen=True)
class ELIspot:
    num_peptides_per_pool: int
    num_coverage: int
    num_processes: int
    peptide_ids: List[str] = field(default_factory=list, repr=False)
    __dummy_peptide_ids: List[str] = field(default_factory=list, repr=False)

    @property
    def num_peptides(self):
        return len(self.peptide_ids)

    @property
    def num_dummy_peptides(self):
        return len(self.__dummy_peptide_ids)

    def __post_init__(self):
        # Step 1. Make sure the number of peptides is bigger than the number of peptides per pool.
        if self.num_peptides < self.num_peptides_per_pool:
            logger.error('Number of peptides per pool is bigger than the total number of peptides.')
            exit(1)

        # Step 2. Add dummy peptides if the number of peptides is not divisible by the number of peptides per pool
        remainder = self.num_peptides % self.num_peptides_per_pool
        if remainder != 0:
            # Get new number of peptides (increased to a divisible number)
            num_dummy_peptides = self.num_peptides_per_pool - remainder
            dummy_peptide_id_idx = 1
            for _ in range(0, num_dummy_peptides):
                while True:
                    dummy_peptide_id = 'dummy_peptide_%i' % dummy_peptide_id_idx
                    if dummy_peptide_id not in self.__dummy_peptide_ids and dummy_peptide_id not in self.peptide_ids:
                        self.__dummy_peptide_ids.append(dummy_peptide_id)
                        break
                    else:
                        dummy_peptide_id_idx += 1

    def generate_configuration(
            self,
            disallowed_peptide_pairs: List[Tuple[str,str]] = []
    ) -> pd.DataFrame:
        """
        Generates an ELIspot assay configuration.

        Parameters
        ----------
        disallowed_peptide_pairs    :   List of tuples (peptide ID, peptide ID).

        Returns
        -------
        status                      :   Status (if optimal, returns cp_model.OPTIMAL).
        df_configuration            :   DataFrame with the following columns:
                                        'pool_id'
                                        'peptide_id'
        """
        # Step 1. Calculate the number of pools per coverage
        num_pools_per_coverage = int((self.num_peptides + self.num_dummy_peptides) / self.num_peptides_per_pool)
        pool_ids = ['pool_' + str(i) for i in range(1, num_pools_per_coverage + 1)]
        coverage_ids = ['coverage_' + str(i) for i in range(1, self.num_coverage + 1)]

        # Step 2. Construct a constraint programming model
        model = cp_model.CpModel()
        data_dict = {
            'coverage_id': [],
            'pool_id': [],
            'peptide_id': [],
            'bool_variable': []
        }
        # Step 3. Initialize the constraint programming model dictionary
        peptide_ids = self.peptide_ids + self.__dummy_peptide_ids
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
                                    (df['peptide_id'] == curr_peptide_id), :]
                model.Add(sum(df_matched['bool_variable'].values.tolist()) == 1)

        # Constraint 2. Each pool in each coverage has exactly the same number of peptides
        for curr_coverage_id in df['coverage_id'].unique():
            for curr_pool_id in df['pool_id'].unique():
                df_matched = df.loc[(df['coverage_id'] == curr_coverage_id) &
                                    (df['pool_id'] == curr_pool_id), :]
                model.Add(sum(df_matched['bool_variable'].values.tolist()) == self.num_peptides_per_pool)

        # Constraint 3. No two peptides are in the same pool more than once
        # At the same time, apply constraints for disallowed peptide pairs
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
            if (peptide_id_1, peptide_id_2) in disallowed_peptide_pairs or \
                    (peptide_id_2, peptide_id_1) in disallowed_peptide_pairs:
                # Pairs cannot appear together in the same pool
                model.Add(sum(peptide_pair_bool_variables) == 0)
            else:
                # All pairs can appear together in the same pool at most once
                model.Add(sum(peptide_pair_bool_variables) <= 1)

        # Step 6. Solve
        logger.info('CP solver started.')
        solver = cp_model.CpSolver()
        solver.parameters.num_search_workers = self.num_processes
        solver.enumerate_all_solutions = False
        status = solver.Solve(model)
        logger.info('CP solver finished.')

        if status == cp_model.OPTIMAL:
            logger.info('Solution is optimal.')

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
                        curr_pool_id = 'pool_' + str(
                            num_pools_per_coverage * (curr_coverage_id_int - 1) + int(curr_pool_id.split('_')[1]))

                    solutions_data['coverage_id'].append(curr_coverage_id)
                    solutions_data['pool_id'].append(curr_pool_id)
                    solutions_data['peptide_id'].append(curr_peptide_id)

            df_configuration = pd.DataFrame(solutions_data)
            df_configuration = df_configuration[df_configuration['peptide_id'].isin(self.__dummy_peptide_ids) == False]
            # df_configuration = df_configuration.drop(columns=['coverage_id'])
        else:
            logger.info('Solution is not optimal.')
            df_configuration = pd.DataFrame()

        return status, df_configuration

    def recycle_configuration(
            self,
            df_template_configuration: pd.DataFrame,
            disallowed_peptide_pairs: List[Tuple[str, str]] = []
    ) -> pd.DataFrame:
        """
        Recycle a pre-computed configuration to scale configuration generation
        for a bigger number of peptides.

        Parameters
        ----------
        df_template_configuration   :   DataFrame of a pre-computed configuration
                                        with the following columns:
                                        'pool_id'
                                        'peptide_id'
        disallowed_peptide_pairs    :   List of tuples (peptide ID, peptide ID).

        Returns
        -------
        df_configuration            :   DataFrame with the following columns:
                                        'pool_id'
                                        'peptide_id'
        """
        # Step 1. Calculate the number of batches
        peptide_ids_template = list(df_template_configuration['peptide_id'].unique())
        num_peptides_template = len(list(df_template_configuration['peptide_id'].unique()))
        num_pools_template = max(int(i.split('_')[1]) for i in list(df_template_configuration['pool_id'].unique()))
        num_batches = math.ceil(self.num_peptides / num_peptides_template)

        # Step 2. Create all necessary combinations of template peptide ID and batch ID
        counter = 0
        peptide_batch_ids = []
        for batch_id in range(1, num_batches + 1):
            for peptide_id in peptide_ids_template:
                peptide_batch_ids.append((peptide_id, batch_id))
                counter += 1
                if counter == self.num_peptides:
                    break
            if counter == self.num_peptides:
                break

        # Step 3. First map the disallowed peptide pairs
        data = {
            'template_peptide_id': [],
            'peptide_id': [],
            'batch_id': []
        }
        lookup_dict = {}
        for pair in disallowed_peptide_pairs:
            peptide_id_1 = pair[0]
            peptide_id_2 = pair[1]
            if peptide_id_1 not in lookup_dict.keys():
                template_peptide_id_1 = peptide_batch_ids[0][0]
                batch_id_1 = peptide_batch_ids[0][1]
                peptide_batch_ids.pop(0)
                lookup_dict[peptide_id_1] = batch_id_1
                data['template_peptide_id'].append(template_peptide_id_1)
                data['peptide_id'].append(peptide_id_1)
                data['batch_id'].append(batch_id_1)
            else:
                batch_id_1 = lookup_dict[peptide_id_1]

            if peptide_id_2 not in lookup_dict.keys():
                for i in peptide_batch_ids:
                    if peptide_batch_ids[i][1] != batch_id_1:
                        break
                template_peptide_id_2 = peptide_batch_ids[i][0]
                batch_id_2 = peptide_batch_ids[i][1]
                peptide_batch_ids.pop(i)
                lookup_dict[peptide_id_2] = batch_id_2
                data['template_peptide_id'].append(template_peptide_id_2)
                data['peptide_id'].append(batch_id_2)
                data['batch_id'].append(batch_id_2)

        # Step 4. Map the peptides in the template configuration to self.peptide_ids
        for peptide_id in self.peptide_ids:
            if peptide_id not in lookup_dict.keys():
                peptide_batch_id = peptide_batch_ids[0]
                peptide_batch_ids.pop(0)
                data['template_peptide_id'].append(peptide_batch_id[0])
                data['peptide_id'].append(peptide_id)
                data['batch_id'].append(peptide_batch_id[1])
        df_peptides = pd.DataFrame(data)

        # Step 5. Recycle pre-computed configuration
        data = {
            'pool_id': [],
            'peptide_id': []
        }
        for curr_batch_id in range(1, num_batches + 1):
            for index, value in df_template_configuration.iterrows():
                curr_template_pool_id = value['pool_id']
                curr_template_peptide_id = value['peptide_id']
                pool_id = 'pool_%i' % (int(curr_template_pool_id.split('_')[1]) + (curr_batch_id - 1) * num_pools_template)
                peptide_id = df_peptides.loc[
                    (df_peptides['template_peptide_id'] == curr_template_peptide_id) &
                    (df_peptides['batch_id'] == curr_batch_id)
                ,'peptide_id'].values.tolist()[0]
                data['pool_id'].append(pool_id)
                data['peptide_id'].append(peptide_id)
        df_configuration = pd.DataFrame(data)
        return df_configuration

    @staticmethod
    def assign_well_ids(
            df_configuration: pd.DataFrame,
            plate_type: str
    ) -> pd.DataFrame:
        """
        Assigns plate and well IDs to an ELIspot configuration.

        Parameters
        ----------
        df_configuration    :   DataFrame with the following columns:
                                'pool_id',
                                'peptide_id'
        plate_type          :   Plate type (allowed values: '96-well plate').

        Returns
        -------
        df_configuration    :   DataFrame with the following columns:
                                'pool_id',
                                'peptide_id'
                                'plate_id'
                                'well_id'
        """
        if plate_type == PlateTypes.PLATE_96_WELLS:
            def get_96_well_plate_ids():
                row_prefixes = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
                col_prefixes = range(1, 13)
                return ['%s%s' % (i[0], i[1]) for i in list(itertools.product(row_prefixes, col_prefixes))]

            curr_plate_id = 1
            curr_well_ids = get_96_well_plate_ids()
            pool_well_ids_dict = {
                'pool_id': [],
                'plate_id': [],
                'well_id': []
            }
            for pool_id in sorted(df_configuration['pool_id'].unique()):
                if len(curr_well_ids) == 0:
                    curr_well_ids = get_96_well_plate_ids()
                    curr_plate_id += 1
                curr_well_id = curr_well_ids[0]
                curr_well_ids.pop(0)
                pool_well_ids_dict['pool_id'].append(pool_id)
                pool_well_ids_dict['plate_id'].append(curr_plate_id)
                pool_well_ids_dict['well_id'].append(curr_well_id)
            df_pool_wells_ids = pd.DataFrame(pool_well_ids_dict)
            df_configuration = pd.merge(df_configuration, df_pool_wells_ids, on='pool_id')
            return df_configuration

    @staticmethod
    def identify_hit_peptides(
            hit_pool_ids: List[str],
            df_configuration: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Identifies hit peptide IDs given read-outs from an ELIspot experiment.

        Parameters
        ----------
        hit_pool_ids                    :   Hit pool IDs.
        df_configuration                :   DataFrame of ELIspot configuration.
                                            Expected columns:
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
        # Step 1. Get the configuration's maximum coverage
        configuration_max_coverage = -1
        for coverage_id in df_configuration['coverage_id'].unique():
            num_coverage = int(coverage_id.split('_')[1])
            if num_coverage > configuration_max_coverage:
                configuration_max_coverage = num_coverage

        # Step 2. Identify hit peptide IDs
        hit_peptides_dict = defaultdict(list)
        for curr_pool_id in hit_pool_ids:
            curr_hit_peptide_ids = df_configuration.loc[df_configuration['pool_id'] == curr_pool_id, 'peptide_id'].values.tolist()
            for curr_hit_peptide_id in curr_hit_peptide_ids:
                hit_peptides_dict[curr_hit_peptide_id].append(curr_pool_id)

        # Step 3. Identify coverage and pool IDs for each hit peptide
        data = {
            'peptide_id': [],
            'pool_ids': [],
            'num_coverage': []
        }
        for key, value in hit_peptides_dict.items():
            data['peptide_id'].append(key)
            data['pool_ids'].append(';'.join(value))
            data['num_coverage'].append(len(value))
        df_hits = pd.DataFrame(data)

        # Step 4. Identify peptide maximum coverage
        hit_peptide_max_coverage = df_hits['num_coverage'].max()
        df_hits_max = df_hits[df_hits['num_coverage'] == hit_peptide_max_coverage]
        if len(df_hits_max) == 0:
            logger.info('Returning as there are no peptides with the desired hit coverage (%ix).' % configuration_max_coverage)
            return pd.DataFrame()

        # Step 5. Identify hit pool IDs and the associated peptide IDs
        hit_pool_ids_dict = defaultdict(list) # key = pool ID, value = list of peptide IDs
        for index, value in df_hits_max.iterrows():
            peptide_id = value['peptide_id']
            pool_ids = value['pool_ids'].split(';')
            for pool_id in pool_ids:
                hit_pool_ids_dict[pool_id].append(peptide_id)

        # Step 6. For the peptides that have the maximum coverage,
        # identify second-round assay peptides
        second_round_assay_peptide_ids = set()
        for peptide_id in df_hits_max['peptide_id'].unique():
            pool_ids = df_configuration.loc[df_configuration['peptide_id'] == peptide_id,'pool_id'].values.tolist()
            unique = False
            for pool_id in pool_ids:
                if len(hit_pool_ids_dict[pool_id]) == 1:
                    unique = True
            if not unique:
                second_round_assay_peptide_ids.add(peptide_id)

        # Step 7. Identify hit peptide IDs.
        data = {
            'peptide_id': [],
            'deconvolution_result': []
        }
        for peptide_id in df_hits_max['peptide_id'].unique():
            if peptide_id not in second_round_assay_peptide_ids:
                data['peptide_id'].append(peptide_id)
                data['deconvolution_result'].append(DeconvolutionResults.HIT)
            else:
                data['peptide_id'].append(peptide_id)
                data['deconvolution_result'].append(DeconvolutionResults.CANDIDATE_HIT)
        df_deconvolution = pd.DataFrame(data)
        df_hits_max = pd.merge(df_hits_max, df_deconvolution, on=['peptide_id'])
        return df_hits_max

    @staticmethod
    def verify_configuration(
            df_configuration: pd.DataFrame,
            num_coverage: int
    ) -> bool:
        """
        Verifies whether a given ELIspot configuration satisfies the following constraints:
        1. Each peptide is in 'num_coverage' number of different pools.
        2. Each peptide is in exactly one unique combination of pool IDs.

        Parameters
        ----------
        df_configuration    :   DataFrame of ELIspot configuration.
                                Expected columns:
                                'coverage_id'
                                'pool_id'
                                'peptide_id'
        num_coverage        :   Coverage.

        Returns
        -------
        is_valid            :   True if the input configuration meets all desired constraints.
                                False otherwise.
        """
        # Step 1. Check if each peptide is in 'num_coverage' number of different pools.
        for peptide_id in list(df_configuration['peptide_id'].unique()):
            pool_ids = (df_configuration.loc[df_configuration['peptide_id'] == peptide_id, 'pool_id'].unique())
            if len(pool_ids) != num_coverage:
                logger.info('%s is in %i different pools (expected: %i).' % (peptide_id, len(pool_ids), num_coverage))
                return False

        # Step 2. Check that each peptide belongs in exactly one unique combination of pool IDs
        pool_id_combinations = set()
        for peptide_id in list(df_configuration['peptide_id'].unique()):
            pool_ids = list(df_configuration.loc[df_configuration['peptide_id'] == peptide_id, 'pool_id'].unique())
            pool_ids = sorted(pool_ids)
            pool_id_combinations.add(','.join(pool_ids))
        if len(pool_id_combinations) != len(df_configuration['peptide_id'].unique()):
            logger.info("%i unique combinations of pool IDs. "
                        "Total number of peptides: %i. "
                        "Each peptide is expected to be assigned a unique combination of pools." %
                        (len(pool_id_combinations), len(df_configuration['peptide_id'].unique())))
            return False

        # Step 3. Configuration neets all constraints
        logger.info("Configuration meets all ACE constraints.")
        return True

