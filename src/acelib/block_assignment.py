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
The purpose of this python3 script is to implement the BlockAssignment dataclass.
"""


import copy
import pandas as pd
import math
import random
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations, product
from typing import Dict, List, Tuple, Type
from .constants import PlateTypes
from .logger import get_logger
from .utilities import convert_peptides_to_dataframe
from .types import Assignments, Peptides, PeptideId, PeptidePairs, PoolId, PlateId, WellId


logger = get_logger(__name__)


@dataclass
class BlockAssignment:
    assignments: Assignments = field(default_factory=dict)
    plate_ids: Dict[PoolId, Tuple[PlateId, WellId]] = field(default_factory=dict)

    @property
    def num_pools(self) -> int:
        df_assignments = self.to_dataframe()
        return len(list(df_assignments['pool_id'].unique()))

    @property
    def peptide_ids(self) -> List[PeptideId]:
        df_assignments = self.to_dataframe()
        return list(df_assignments['peptide_id'].unique())

    @property
    def pooled_peptide_pairs(self) -> PeptidePairs:
        """
        Returns a list of peptide pairs that are pooled together.

        Returns
        -------
        peptide_pairs   :   List of all peptide pairs that appear together in a pool.
        """
        peptide_pairs = []
        df_assignments = self.to_dataframe()
        for name, group in df_assignments.groupby('pool_id'):
            for peptide_id_1, peptide_id_2 in combinations(group['peptide_id'], r=2):
                peptide_pairs.append((peptide_id_1, peptide_id_2))
        return peptide_pairs

    def add_peptide(
            self,
            coverage: int,
            pool: int,
            peptide_id: str,
            peptide_sequence: str,
    ):
        """
        Add a peptide.

        Parameters
        ----------
        coverage            :   Coverage ID.
        pool                :   Pool ID.
        peptide_id          :   Peptide ID.
        peptide_sequence    :   Peptide sequence.
        """
        if coverage not in self.assignments.keys():
            self.assignments[coverage] = {}
        if pool not in self.assignments[coverage].keys():
            self.assignments[coverage][pool] = []
        self.assignments[coverage][pool].append((peptide_id, peptide_sequence))

    def count_violations(self) -> int:
        """
        Counts the number of violations (i.e. number of peptides with
        non-unique pool assignment).

        Returns
        -------
        num_violations      :   Number of violations.
        """
        df_assignments = self.to_dataframe()
        pool_ids_peptides_dict = defaultdict(list)
        for peptide_id in list(df_assignments['peptide_id'].unique()):
            pool_ids = list(df_assignments.loc[df_assignments['peptide_id'] == peptide_id, 'pool_id'].unique())
            pool_ids = sorted(pool_ids)
            pool_ids_peptides_dict[','.join([str(i) for i in pool_ids])].append(peptide_id)

        num_violations = 0
        for key, value in pool_ids_peptides_dict.items():
            if len(value) > 1:
                for peptide_id in value:
                    num_violations += 1
        return num_violations

    def to_dataframe(self) -> pd.DataFrame:
        data = {
            'coverage_id': [],
            'pool_id': [],
            'peptide_id': [],
            'peptide_sequence': [],
            'plate_id': [],
            'well_id': []
        }
        for coverage_id in self.assignments.keys():
            for pool_id in self.assignments[coverage_id].keys():
                for peptide_id, peptide_sequence in self.assignments[coverage_id][pool_id]:
                    data['coverage_id'].append(coverage_id)
                    data['pool_id'].append(pool_id)
                    data['peptide_id'].append(peptide_id)
                    data['peptide_sequence'].append(peptide_sequence)
                    if pool_id in self.plate_ids.keys():
                        data['plate_id'].append(self.plate_ids[pool_id][0])
                        data['well_id'].append(self.plate_ids[pool_id][1])
                    else:
                        data['plate_id'].append('')
                        data['well_id'].append('')
        return pd.DataFrame(data)

    def is_optimal(
            self,
            num_coverage: int,
            num_peptides_per_pool: int,
            verbose: bool = True
    ) -> bool:
        """
        Verifies whether a given ELISpot assignment satisfies the following constraints:
        1. Each peptide is in 'num_coverage' number of different pools.
        2. Each peptide is in exactly one unique combination of pool IDs.
        3. There is an optimal (minimal) number of pools.

        Parameters
        ---------
        num_coverage            :   Coverage.
        num_peptides_per_pool   :   Number of peptides per pool.
        verbose                 :   Verbose.

        Returns
        -------
        is_optimal              :   True if the input configuration meets all desired criteria.
                                    False otherwise.
        """
        df_assignments = self.to_dataframe()

        # Step 1. Check if each peptide is in 'num_coverage' number of different pools.
        constraint_1_bool = True
        for peptide_id in list(df_assignments['peptide_id'].unique()):
            pool_ids = (df_assignments.loc[df_assignments['peptide_id'] == peptide_id, 'pool_id'].unique())
            if len(pool_ids) != num_coverage:
                if verbose:
                    logger.info('Assignment does not meet constraint #1: peptide %s is in %i different pools (expected: %i).' %
                                (peptide_id, len(pool_ids), num_coverage))
                constraint_1_bool = False
        if constraint_1_bool:
            if verbose:
                logger.info('Assignment meets constraint #1: each peptide is in %i different pools.' % num_coverage)

        # Step 2. Check that each peptide belongs to exactly one unique combination of pool IDs
        constraint_2_bool = True
        pool_ids_peptides_dict = defaultdict(list)
        for peptide_id in list(df_assignments['peptide_id'].unique()):
            pool_ids = list(df_assignments.loc[df_assignments['peptide_id'] == peptide_id, 'pool_id'].unique())
            pool_ids = sorted(pool_ids)
            pool_ids_peptides_dict[','.join([str(i) for i in pool_ids])].append(peptide_id)
        for key, value in pool_ids_peptides_dict.items():
            if len(value) > 1:
                if verbose:
                    logger.info("Assignment does not meet constraint #2. Pools %s have the following peptides:" % key)
                    for peptide_id in value:
                        logger.info(peptide_id)
                constraint_2_bool = False
        if constraint_2_bool:
            if verbose:
                logger.info('Assignment meets constraint #2: each peptide belongs to exactly one unique combination of pool IDs.')

        # Step 3. Check that there is an optimal number of pools
        constraint_3_bool = True
        num_pools = math.ceil(len(df_assignments['peptide_id'].unique()) / num_peptides_per_pool) * num_coverage
        if len(df_assignments['pool_id'].unique()) != num_pools:
            num_extra_pools = len(df_assignments['pool_id'].unique()) - num_pools
            if verbose:
                logger.info('Assignment does not meet constraint #3: %i extra pool(s) than the minimum possible number of pools (%i).' %
                            (num_extra_pools, num_pools))
            constraint_3_bool = False
        if constraint_3_bool:
            if verbose:
                logger.info('Assignment meets constraint #3: there is an optimal (minimal) number of pools (%i).' % num_pools)

        return constraint_1_bool & constraint_2_bool & constraint_3_bool

    def shuffle_pool_ids(self):
        """
        Shuffles

        Returns
        -------

        """
        # Step 1. Shuffle pool IDs
        curr_pool_ids = list(self.to_dataframe()['pool_id'].unique())
        new_pool_ids = list(self.to_dataframe()['pool_id'].unique())
        random.shuffle(new_pool_ids)

        # Step 2. Create a dictionary of old and new pool IDs
        new_pool_ids_dict = {}
        for i in range(0, len(curr_pool_ids)):
            new_pool_ids_dict[curr_pool_ids[i]] = new_pool_ids[i]

        # Step 3. Reassign pool IDs
        old_assignments = self.assignments
        self.assignments = {}
        for coverage, pool_dict in old_assignments.items():
            for pool, peptide_ids in pool_dict.items():
                for peptide_id, peptide_sequence in peptide_ids:
                    self.add_peptide(
                        coverage=coverage,
                        pool=new_pool_ids_dict[pool],
                        peptide_id=peptide_id,
                        peptide_sequence=peptide_sequence
                    )

    @staticmethod
    def assign_well_ids(
            df_assignment: pd.DataFrame,
            plate_type: str
    ) -> pd.DataFrame:
        """
        Assigns plate and well IDs to an ELISpot configuration.

        Parameters
        ----------
        df_assignment       :   DataFrame with the following columns:
                                'coverage_id'
                                'pool_id',
                                'peptide_id'
                                'peptide_sequence'
        plate_type          :   Plate type (allowed values: '96-well plate').

        Returns
        -------
        df_assignment       :   DataFrame with the following columns:
                                'coverage_id'
                                'pool_id',
                                'peptide_id'
                                'peptide_sequence'
                                'plate_id'
                                'well_id'
        """
        if plate_type == PlateTypes.PLATE_96_WELLS:
            def get_96_well_plate_ids():
                row_prefixes = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
                col_prefixes = range(1, 13)
                return ['%s%s' % (i[0], i[1]) for i in list(product(row_prefixes, col_prefixes))]

            curr_plate_id = 1
            curr_well_ids = get_96_well_plate_ids()
            pool_well_ids_dict = {
                'pool_id': [],
                'plate_id': [],
                'well_id': []
            }
            for pool_id in sorted(df_assignment['pool_id'].unique()):
                if len(curr_well_ids) == 0:
                    curr_well_ids = get_96_well_plate_ids()
                    curr_plate_id += 1
                curr_well_id = curr_well_ids[0]
                curr_well_ids.pop(0)
                pool_well_ids_dict['pool_id'].append(pool_id)
                pool_well_ids_dict['plate_id'].append(curr_plate_id)
                pool_well_ids_dict['well_id'].append(curr_well_id)
            df_pool_wells_ids = pd.DataFrame(pool_well_ids_dict)
            df_assignment.drop(columns=['plate_id', 'well_id'], inplace=True)
            df_configuration = pd.merge(df_assignment, df_pool_wells_ids, on='pool_id')
            return df_configuration
        else:
            logger.error('Unsupported plate_type: %s' % plate_type)
            exit(1)

    @staticmethod
    def generate_single_coverage_block_assignment(
            peptides: Peptides,
            preferred_peptide_pairs: PeptidePairs,
            num_peptides_per_pool: int,
            coverage: int = 1
    ) -> 'BlockAssignment':
        """
        Generate a one-coverage block assignment given a list of preferred peptide pairs.

        Parameters
        ----------
        peptides                :   Peptides.
        preferred_peptide_pairs :   PeptidePairs.
        num_peptides_per_pool   :   Number of peptides per pool.
        coverage                :   Coverage (default: 1).

        Returns
        -------
        block_assignment        :   BlockAssignment object.
        """
        df_peptides = convert_peptides_to_dataframe(peptides=peptides)

        # Step 1. Initialize pools
        pools = {}
        num_pools = math.ceil(len(df_peptides['peptide_id'].unique()) / num_peptides_per_pool)
        for i in range(1, num_pools + 1):
            pools[i] = []

        # Step 2. Compute transitive neighbors
        peptide_neighbors = compute_transitive_neighbors(peptide_pairs=preferred_peptide_pairs)

        # Step 3. Assign pool IDs for peptide neighbors
        preferred_peptide_ids = []
        for peptide_neighbor in peptide_neighbors:
            for peptide_id in peptide_neighbor:
                preferred_peptide_ids.append(peptide_id)
                for pool_id in pools.keys():
                    if len(pools[pool_id]) < num_peptides_per_pool:
                        pools[pool_id].append(peptide_id)
                        break

        # Step 4. Assign pool IDs for remaining peptide IDs
        for peptide_id in df_peptides['peptide_id'].unique():
            if peptide_id not in preferred_peptide_ids:
                for curr_pool_idx in range(1, num_pools + 1):
                    if len(pools[curr_pool_idx]) < num_peptides_per_pool:
                        pools[curr_pool_idx].append(peptide_id)
                        break

        # Step 5. Assign pool IDs for all peptides
        block_assignment = BlockAssignment()
        for key, val in pools.items():
            for peptide_id in val:
                peptide_sequence = df_peptides.loc[df_peptides['peptide_id'] == peptide_id, 'peptide_sequence'].values[
                    0]
                block_assignment.add_peptide(
                    coverage=coverage,
                    pool=key,
                    peptide_id=peptide_id,
                    peptide_sequence=peptide_sequence
                )
        return block_assignment

    @staticmethod
    def load_from_dataframe(df_assignments: pd.DataFrame) -> 'BlockAssignment':
        """
        Loads assignments from a pd.DataFrame and returns a BlockAssignment object.

        Parameters
        ----------
        df_assignments      :   pd.DataFrame with the following columns:
                                'coverage_id'
                                'pool_id'
                                'peptide_id'
                                'peptide_sequence'

        Returns
        -------
        block_assignment    :   BlockAssignment object.
        """
        block_assignment = BlockAssignment()
        for index, row in df_assignments.iterrows():
            coverage_id = int(row['coverage_id'])
            pool_id = int(row['pool_id'])
            peptide_id = str(row['peptide_id'])
            peptide_sequence = str(row['peptide_sequence'])
            block_assignment.add_peptide(
                coverage=coverage_id,
                pool=pool_id,
                peptide_id=peptide_id,
                peptide_sequence=peptide_sequence
            )
            if 'plate_id' in row and 'well_id' in row:
                plate_id = row['plate_id']
                well_id = row['well_id']
                block_assignment.plate_ids[row['pool_id']] = (plate_id, well_id)
        return block_assignment

    @staticmethod
    def load_golfy_assignment(
            golfy_assignment: Dict,
            df_peptides: pd.DataFrame
    ) -> 'BlockAssignment':
        """
        Converts golfy results into a DataFrame.

        Parameters
        ----------
        golfy_assignment    :   Dictionary.
        df_peptides         :   pd.DatFrame with the following columns:
                                'peptide_id'
                                'peptide_index'
                                'peptide_sequence'

        Returns
        -------
        block_assignment    :   BlockAssignment object.
        """
        block_assignment = BlockAssignment()
        curr_pool_id = 1
        for key, value in golfy_assignment.items():
            curr_coverage_id = key + 1
            for key2, value2 in value.items():
                for p in value2:
                    df_curr_peptide = df_peptides.loc[df_peptides['peptide_index'] == p, :]
                    curr_peptide_id = df_curr_peptide['peptide_id'].values[0]
                    curr_peptide_sequence = df_curr_peptide['peptide_sequence'].values[0]
                    block_assignment.add_peptide(
                        coverage=curr_coverage_id,
                        pool=curr_pool_id,
                        peptide_id=curr_peptide_id,
                        peptide_sequence=curr_peptide_sequence
                    )
                curr_pool_id += 1
        return block_assignment

    @staticmethod
    def read_excel_file(
            excel_file: str,
            sheet_name: str = 'block_assignment'
    ) -> 'BlockAssignment':
        """
        Reads an Excel file and returns a BlockAssignment object.

        Parameters
        ----------
        excel_file          :   Excel file.
        sheet_name          :   Sheet name (default: 'block_assignment').

        Returns
        -------
        block_assignment    :   BlockAssignment object.
        """
        df = pd.read_excel(excel_file, sheet_name=sheet_name)
        block_assignment = BlockAssignment()
        for index, row in df.iterrows():
            block_assignment.add_peptide(
                coverage=row['coverage_id'],
                pool=row['pool_id'],
                peptide_id=row['peptide_id'],
                peptide_sequence=row['peptide_sequence']
            )
            if 'plate_id' in row and 'well_id' in row:
                plate_id = row['plate_id']
                well_id = row['well_id']
                block_assignment.plate_ids[row['pool_id']] = (plate_id, well_id)
        return block_assignment

    @staticmethod
    def update_ids(
            assignments: Assignments,
            start_pool_num: int,
            start_coverage_num: int
    ) -> Assignments:
        """
        Updates IDs of an ELISpot configuration.

        Parameters
        ----------
        assignments             :   Assignments.
        start_pool_num          :   Start pool number (e.g. 1).
        start_coverage_num      :   Start coverage number (e.g. 1).

        Returns
        -------
        new_assignments         :   Assignments.
        """
        # Step 1. Get all coverages and pools
        old_coverages = sorted(assignments.keys())
        old_pools = []
        for coverage in old_coverages:
            for pool in assignments[coverage].keys():
                old_pools.append(pool)
        old_pools = sorted(old_pools)

        # Step 2. Create a dictionary of old and new coverages
        coverages_dict = {} # key = old coverage, value = new coverage
        curr_coverage_num = start_coverage_num
        for coverage in old_coverages:
            coverages_dict[coverage] = curr_coverage_num
            curr_coverage_num += 1

        # Step 3. Create a dictionary of old and new pools
        pools_dict = {} # key = old pool, value = new pool
        curr_pool_num = start_pool_num
        for pool in old_pools:
            pools_dict[pool] = curr_pool_num
            curr_pool_num += 1

        # Step 4. Prepare new assignments
        new_assignments = {}
        for coverage in assignments.keys():
            new_assignments[coverages_dict[coverage]] = {}
            for pool in assignments[coverage].keys():
                new_assignments[coverages_dict[coverage]][pools_dict[pool]] = []

        # Step 5. Reassign peptides into new coverages and new pools
        for coverage in assignments.keys():
            new_assignments[coverages_dict[coverage]] = {}
            for pool in assignments[coverage].keys():
                new_assignments[coverages_dict[coverage]][pools_dict[pool]] = []
                for peptide in assignments[coverage][pool]:
                    new_assignments[coverages_dict[coverage]][pools_dict[pool]].append(peptide)

        return new_assignments

    @staticmethod
    def merge(block_assignments: List['BlockAssignment']) -> 'BlockAssignment':
        """
        Merges a list of

        Parameters
        ----------
        block_assignments   :   List of BlockAssignment objects.

        Returns
        -------
        block_assignment    :   BlockAssignment object.
        """
        block_assignment = BlockAssignment()
        for block_assignment_ in block_assignments:
            df_assignment = block_assignment_.to_dataframe()
            for index, row in df_assignment.iterrows():
                block_assignment.add_peptide(
                    coverage=row['coverage_id'],
                    pool=row['pool_id'],
                    peptide_id=row['peptide_id'],
                    peptide_sequence=row['peptide_sequence']
                )
        return block_assignment

    @staticmethod
    def minimize_violations(
            block_assignments: List['BlockAssignment'],
            shuffle_iters: int,
            verbose: bool = True
    ) -> List['BlockAssignment']:
        """
        Minimizes violations (i.e. non-unique pool assignment) in a list of
        block assignments by shuffling pool IDs.

        Parameters
        ----------
        block_assignments   :   List of BlockAssignment objects.
        shuffle_iters       :   Number of iterations to shuffle pool IDs.

        Returns
        -------
        block_assignments   :   List of BlockAssignment objects.
        """
        min_violations = BlockAssignment.merge(block_assignments=block_assignments).count_violations()
        curr_block_assignments = copy.deepcopy(block_assignments)
        best_block_assignments = copy.deepcopy(block_assignments)
        for _ in range(0, shuffle_iters):
            random_idx = random.choice(list(range(0, len(curr_block_assignments))))
            curr_block_assignments[random_idx].shuffle_pool_ids()
            curr_num_violations = BlockAssignment.merge(block_assignments=curr_block_assignments).count_violations()
            if curr_num_violations < min_violations:
                if verbose:
                    logger.info('Found a better assignment: current number of violations: %i, new number of violations: %i' %
                                (min_violations, curr_num_violations))
                best_block_assignments = copy.deepcopy(curr_block_assignments)
                min_violations = curr_num_violations
        return best_block_assignments


def compute_transitive_neighbors(
        peptide_pairs: PeptidePairs
) -> List[List[PeptideId]]:
    """
    Computes transitive peptide neighbors.

    Parameters
    ----------
    peptide_pairs       :   List of tuples (peptide ID, peptide ID).

    Returns
    -------
    peptide_neighbors   :   List of list of peptide IDs.
    """
    # Step 1. Create a dictionary to store peptide relationships
    peptide_dict = {}
    for peptide_id_1, peptide_id_2 in peptide_pairs:
        if peptide_id_1 not in peptide_dict:
            peptide_dict[peptide_id_1] = set()
        if peptide_id_2 not in peptide_dict:
            peptide_dict[peptide_id_2] = set()
        peptide_dict[peptide_id_1].add(peptide_id_2)
        peptide_dict[peptide_id_2].add(peptide_id_1)

    # Step 2. Find transitive peptide IDs
    transitive_peptides = []
    for peptide_id in peptide_dict:
        stack = list(peptide_dict[peptide_id])
        visited = set()
        while stack:
            curr_peptide_id = stack.pop()
            if curr_peptide_id not in visited:
                visited.add(curr_peptide_id)
                if curr_peptide_id in peptide_dict:
                    stack.extend(peptide_dict[curr_peptide_id])
        transitive_peptides.append(list(visited))

    # Step 3. De-duplicate peptide IDs
    transitive_peptides_ = []
    for i in transitive_peptides:
        duplicate = False
        for j in transitive_peptides_:
            if set(i) == j:
                duplicate = True
                break
        if not duplicate:
            transitive_peptides_.append(set(i))

    return [list(p) for p in transitive_peptides_]


