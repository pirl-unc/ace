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
The purpose of this python3 script is to implement the BlockDesign dataclass.
"""


import math
import pandas as pd
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations
from ortools.sat.python import cp_model
from typing import Dict, List, Type
from .block_assignment import BlockAssignment
from .types import *
from .logger import get_logger


logger = get_logger(__name__)


@dataclass
class BlockDesign:
    peptides: Peptides
    num_peptides_per_pool: int
    num_coverage: int
    max_peptides_per_block: int
    disallowed_peptide_pairs: PeptidePairs = field(default_factory=list)
    preferred_peptide_pairs: PeptidePairs = field(default_factory=list)
    __dummy_peptide_ids: List[PeptideId] = field(default_factory=list, repr=False)
    __peptides_dict: Dict[str, str] = field(default_factory=dict)

    @property
    def all_peptide_ids(self) -> List[PeptideId]:
        return self.peptide_ids + self.__dummy_peptide_ids

    @property
    def num_dummy_peptides(self) -> int:
        return len(self.__dummy_peptide_ids)

    @property
    def num_peptides(self) -> int:
        return len(self.peptides)

    @property
    def num_total_peptides(self) -> int:
        return self.num_peptides + self.num_dummy_peptides

    @property
    def peptide_ids(self) -> List[PeptideId]:
        peptide_ids = []
        for peptide_id, peptide_sequence in self.peptides:
            peptide_ids.append(peptide_id)
        return peptide_ids

    def __post_init__(self):
        # Step 1. Make sure the number of peptides is bigger than the number of peptides per pool.
        if self.num_peptides < self.num_peptides_per_pool:
            logger.error('Number of peptides per pool is bigger than the total number of peptides.')
            exit(1)

        # Step 2. Add dummy peptide IDs
        num_dummy_peptides = self.max_peptides_per_block - self.num_peptides
        dummy_peptide_id_idx = 0
        for _ in range(0, num_dummy_peptides):
            while True:
                dummy_peptide_id_idx += 1
                dummy_peptide_id = 'dummy_peptide_%i' % dummy_peptide_id_idx
                if dummy_peptide_id not in self.peptide_ids:
                    self.__dummy_peptide_ids.append(dummy_peptide_id)
                    break

        # Step 3. Populate peptide dictionary
        for peptide in self.peptides:
            self.__peptides_dict[peptide[0]] = peptide[1]

    def get_peptide_sequence(self, peptide_id: PeptideId) -> PeptideSequence:
        return self.__peptides_dict[peptide_id]

    def generate(
            self,
            random_seed: int,
            num_processes: int,
            verbose: bool = True
    ) -> BlockAssignment:
        """
        Generates ELISpot block design assignments.

        Parameters
        ----------
        random_seed         :   Random seed.
        num_processes       :   Number of processes.
        verbose             :   If True, prints messages.

        Returns
        -------
        block_assignment    :   BlockAssignment object.
        """
        # Step 1. Calculate the number of pools per coverage
        num_pools_per_coverage = int(self.num_total_peptides / self.num_peptides_per_pool)
        pool_ids = list(range(0, num_pools_per_coverage))
        coverage_ids = list(range(0, self.num_coverage))

        # Step 2. Construct a constraint programming model
        model = cp_model.CpModel()
        data_dict = {
            'coverage_id': [],
            'pool_id': [],
            'peptide_id': [],
            'bool_variable': []
        }

        # Step 3. Initialize the constraint programming model dictionary
        var_dict = {}
        for curr_coverage_id in coverage_ids:
            for curr_pool_id in pool_ids:
                for curr_peptide_id in self.all_peptide_ids:
                    curr_bool_var = model.NewBoolVar("%i/%i/%s" % (curr_coverage_id, curr_pool_id, curr_peptide_id))
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
        for peptide_id_1, peptide_id_2 in combinations(self.all_peptide_ids, r=2):
            peptide_pair_bool_variables = []
            for curr_coverage_id in coverage_ids:
                for curr_pool_id in pool_ids:
                    pair_bool_variable = model.NewBoolVar("%i/%s/%s" % (curr_coverage_id, peptide_id_1, peptide_id_2))
                    peptide_1_bool_variable = var_dict[(curr_coverage_id, curr_pool_id, peptide_id_1)]
                    peptide_2_bool_variable = var_dict[(curr_coverage_id, curr_pool_id, peptide_id_2)]

                    # pair_bool_variable has to be 1 if peptide 1 and peptide 2 are paired together
                    model.Add((peptide_1_bool_variable + peptide_2_bool_variable - pair_bool_variable) <= 1)
                    peptide_pair_bool_variables.append(pair_bool_variable)

            # Include boolean variable to apply disallowed peptide pairs or enforced peptide pairs
            if (peptide_id_1, peptide_id_2) in self.disallowed_peptide_pairs or \
                    (peptide_id_2, peptide_id_1) in self.disallowed_peptide_pairs:
                # Pairs cannot appear together in the same pool
                model.Add(sum(peptide_pair_bool_variables) == 0)
            else:
                # All pairs can appear together in the same pool at most once
                model.Add(sum(peptide_pair_bool_variables) <= 1)

        # Step 4. Solve
        solver = cp_model.CpSolver()
        solver.parameters.num_search_workers = num_processes
        solver.enumerate_all_solutions = False
        solver.parameters.random_seed = random_seed
        status = solver.Solve(model)

        if status == cp_model.OPTIMAL:
            if verbose:
                logger.info('\tAn optimal feasible solution was found.')
        elif status == cp_model.FEASIBLE:
            if verbose:
                logger.info('\tA feasible solution was found, but we do not know if it is optimal.')
        elif status == cp_model.INFEASIBLE:
            logger.error('\tThe problem was proven infeasible.')
            exit(1)
        elif status == cp_model.MODEL_INVALID:
            logger.error('\tThe given CpModelProto did not pass the validation step.')
            exit(1)
        elif status == cp_model.UNKNOWN:
            logger.error('\tThe status of the model is unknown because no solution was '
                         'found before something caused the solver to stop, such as a time limit or a memory limit.')
            exit(1)

        # Step 5. Parse solution
        block_assignment = BlockAssignment()
        for curr_bool_variable in df['bool_variable'].values.tolist():
            if solver.Value(curr_bool_variable) == 1:
                curr_bool_variable_elements = str(curr_bool_variable).split("/")
                curr_coverage_id = int(curr_bool_variable_elements[0])
                curr_pool_id = int(curr_bool_variable_elements[1])
                curr_peptide_id = str(curr_bool_variable_elements[2])

                if curr_peptide_id not in self.__dummy_peptide_ids:
                    curr_pool_id = (num_pools_per_coverage * curr_coverage_id) + curr_pool_id
                    curr_peptide_sequence = self.get_peptide_sequence(peptide_id=curr_peptide_id)
                    block_assignment.add_peptide(
                        coverage=curr_coverage_id,
                        pool=curr_pool_id,
                        peptide_id=curr_peptide_id,
                        peptide_sequence=curr_peptide_sequence
                    )

        return block_assignment

    def to_dataframe(self) -> pd.DataFrame:
        peptides = []
        for peptide_id, peptide_sequence in self.peptides:
            peptides.append('%s|%s' % (peptide_id, peptide_sequence))
        disallowed_peptide_pairs = []
        for peptide_id_1, peptide_id_2, in self.disallowed_peptide_pairs:
            disallowed_peptide_pairs.append('%s|%s' % (peptide_id_1, peptide_id_2))
        preferred_peptide_pairs = []
        for peptide_id_1, peptide_id_2, in self.preferred_peptide_pairs:
            preferred_peptide_pairs.append('%s|%s' % (peptide_id_1, peptide_id_2))
        data = {
            'peptides': [';'.join(peptides)],
            'num_peptides_per_pool': [self.num_peptides_per_pool],
            'num_coverage': [self.num_coverage],
            'max_peptides_per_block': [self.max_peptides_per_block],
            'disallowed_peptide_pairs': [';'.join(disallowed_peptide_pairs)],
            'preferred_peptide_pairs': [';'.join(preferred_peptide_pairs)]
        }
        return pd.DataFrame(data)

    @staticmethod
    def compute_num_total_pools(
            num_peptides: int,
            num_peptides_per_design: int,
            num_peptides_per_pool: int,
            num_coverage: int
    ) -> int:
        """
        Computes the total number of pools (theoretical minimum).

        Parameters
        ----------
        num_peptides                :   Number of peptides.
        num_peptides_per_design     :   Number of peptides per design.
        num_peptides_per_pool       :   Number of peptides per pool.
        num_coverage                :   Coverage.

        Returns
        -------
        num_total_pools             :   Total number of pools.
        """
        num_total_pools = 0
        while True:
            num_total_pools += math.ceil(num_peptides_per_design / num_peptides_per_pool) * num_coverage
            num_peptides -= num_peptides_per_design
            if num_peptides <= 0:
                break
        return num_total_pools

    @staticmethod
    def divide_block_design(
            block_design: 'BlockDesign',
            max_peptides_per_block: int,
            max_peptides_per_pool: int,
            verbose: bool = True
    ) -> List[List['BlockDesign']]:
        """
        Divide a block design into multiple computationally tractable designs.

        Parameters
        ----------
        block_design            :   BlockDesign object.
        max_peptides_per_block  :   Maximum number of peptides per block.
        max_peptides_per_pool   :   Maximum number of peptides per pool.
        verbose                 :   If True, prints messages.

        Returns
        -------
        block_designs           :   List of BlockDesign object lists
                                    (block designs in the same list should be
                                     merged into one).
        """
        # Step 1. Divide the number of peptides per pool by two until each
        # one is less than or equal to max_peptides_per_pool
        list_num_peptides_per_pool = [block_design.num_peptides_per_pool]
        while True:
            list_num_peptides_per_pool_ = []
            for i in list_num_peptides_per_pool:
                if i > max_peptides_per_pool:
                    val1 = math.ceil(i / 2)
                    val2 = i - val1
                    list_num_peptides_per_pool_.append(val1)
                    list_num_peptides_per_pool_.append(val2)
                else:
                    list_num_peptides_per_pool_.append(i)

            # Check if every value is less than or equal to max_peptides_per_pool
            list_num_peptides_per_pool = list_num_peptides_per_pool_
            is_done = True
            for i in list_num_peptides_per_pool:
                if i > max_peptides_per_pool:
                    is_done = False
            if is_done:
                break

        # Step 2. Divide the number of peptides proportionately to the
        # number of peptides per pool
        list_num_peptides = []
        num_remaining_peptides = block_design.num_peptides
        for num_peptides_per_pool in list_num_peptides_per_pool:
            num_peptides = math.ceil((num_peptides_per_pool / block_design.num_peptides_per_pool) * block_design.num_peptides)
            if num_remaining_peptides >= num_peptides:
                list_num_peptides.append(num_peptides)
                num_remaining_peptides -= num_peptides
            else:
                list_num_peptides.append(num_remaining_peptides)
                num_remaining_peptides -= num_remaining_peptides

        if verbose:
            logger.info('Divided block design into:')
            for i in range(0, len(list_num_peptides)):
                logger.info('\t%i peptides; %i peptides per pool' %
                            (list_num_peptides[i],
                             list_num_peptides_per_pool[i]))

        # Step 3. Divide the block design into smaller designs such that
        # the number of peptides per design is equal to or smaller than
        # the maximum number of peptides per design
        block_designs = []
        start_peptide_idx = 0
        for i in range(0, len(list_num_peptides_per_pool)):
            num_peptides_per_pool = list_num_peptides_per_pool[i]
            num_peptides = list_num_peptides[i]
            if verbose:
                logger.info('Computing the optimal number of peptides per design (%i peptides; %i peptides per pool).' %
                            (num_peptides, num_peptides_per_pool))
            block_designs_ = []
            if num_peptides > max_peptides_per_block:
                # Figure out the optimal number of peptides per block
                min_peptides = num_peptides_per_pool * num_peptides_per_pool
                max_peptides = max_peptides_per_block
                optimal_num_peptides_per_design = -1
                min_total_pools = -1
                for curr_num_peptides in range(min_peptides, max_peptides + 1):
                    num_total_pools = BlockDesign.compute_num_total_pools(
                        num_peptides=num_peptides,
                        num_peptides_per_design=curr_num_peptides,
                        num_peptides_per_pool=num_peptides_per_pool,
                        num_coverage=block_design.num_coverage
                    )
                    if min_total_pools == -1:
                        min_total_pools = num_total_pools
                        optimal_num_peptides_per_design = curr_num_peptides
                    else:
                        if num_total_pools < min_total_pools:
                            min_total_pools = num_total_pools
                            optimal_num_peptides_per_design = curr_num_peptides
                if verbose:
                    logger.info('\tOptimal number of peptides: %i' % optimal_num_peptides_per_design)
            else:
                optimal_num_peptides_per_design = num_peptides

            last_peptide_idx = start_peptide_idx + num_peptides - 1
            while True:
                end_peptide_idx = start_peptide_idx + optimal_num_peptides_per_design - 1
                if end_peptide_idx > last_peptide_idx:
                    end_peptide_idx = last_peptide_idx
                peptides = block_design.peptides[start_peptide_idx:end_peptide_idx + 1]
                start_peptide_idx = end_peptide_idx + 1
                if verbose:
                    logger.info('\tAppending block design for %i peptides, %i peptides per pool' %
                                (len(peptides), num_peptides_per_pool))
                block_design_ = BlockDesign(
                    peptides=peptides,
                    num_peptides_per_pool=num_peptides_per_pool,
                    num_coverage=block_design.num_coverage,
                    max_peptides_per_block=optimal_num_peptides_per_design,
                    disallowed_peptide_pairs=block_design.disallowed_peptide_pairs,
                    preferred_peptide_pairs=block_design.preferred_peptide_pairs
                )
                block_designs_.append(block_design_)
                if end_peptide_idx == last_peptide_idx:
                    break

            block_designs.append(block_designs_)
        return block_designs

    @staticmethod
    def read_excel_file(
            excel_file: str,
            sheet_name: str = 'block_design'
    ) -> 'BlockDesign':
        """
        Reads an Excel file and returns a BlockDesign object.

        Parameters
        ----------
        excel_file      :   Excel file.
        sheet_name      :   Sheet name.

        Returns
        -------
        block_design    :   BlockDesign object.
        """
        df = pd.read_excel(excel_file, sheet_name=sheet_name)
        peptides = []
        disallowed_peptide_pairs = []
        preferred_peptide_pairs = []
        for peptide in str(df['peptides'].values[0]).split(';'):
            peptides.append((peptide.split('|')[0], peptide.split('|')[1]))
        if not pd.isna(df['disallowed_peptide_pairs'].values[0]):
            peptide_pairs = str(df['disallowed_peptide_pairs'].values[0]).split(';')
            for peptide_pair in peptide_pairs:
                if peptide_pair != '':
                    peptide_id_1 = peptide_pair.split('|')[0]
                    peptide_id_2 = peptide_pair.split('|')[1]
                    disallowed_peptide_pairs.append((peptide_id_1, peptide_id_2))
        if not pd.isna(df['preferred_peptide_pairs'].values[0]):
            peptide_pairs = str(df['preferred_peptide_pairs'].values[0]).split(';')
            for peptide_pair in peptide_pairs:
                if peptide_pair != '':
                    peptide_id_1 = peptide_pair.split('|')[0]
                    peptide_id_2 = peptide_pair.split('|')[1]
                    preferred_peptide_pairs.append((peptide_id_1, peptide_id_2))
        num_peptides_per_pool = df['num_peptides_per_pool'].values[0]
        num_coverage = df['num_coverage'].values[0]
        max_peptides_per_block = df['max_peptides_per_block'].values[0]
        block_design = BlockDesign(
            peptides=peptides,
            num_peptides_per_pool=num_peptides_per_pool,
            num_coverage=num_coverage,
            max_peptides_per_block=max_peptides_per_block,
            disallowed_peptide_pairs=disallowed_peptide_pairs,
            preferred_peptide_pairs=preferred_peptide_pairs
        )
        return block_design

