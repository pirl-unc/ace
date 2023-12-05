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
from golfy import Design
from itertools import combinations, product
from typing import Dict, List, Tuple, Type
from .constants import *
from .logger import get_logger
from .utilities import convert_peptides_to_dataframe
from .types import *


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

    def assign_well_ids(self, plate_size: int):
        """
        Assigns plate and well IDs to an ELISpot configuration.

        Parameters
        ----------
        plate_size     :   Number of wells on plate (allowed values: 24, 48, 96, 384).
        """
        # Step 1. Clear the current self.plate_ids
        self.plate_ids = {} # key = pool ID, value = (plate ID, well ID)

        # Step 2. Prepare well IDs
        def get_24_well_ids():
            row_prefixes = ['A', 'B', 'C', 'D']
            col_prefixes = range(1, 7)
            return ['%s%s' % (i[0], i[1]) for i in list(product(row_prefixes, col_prefixes))]
        def get_48_well_ids():
            row_prefixes = ['A', 'B', 'C', 'D', 'E', 'F']
            col_prefixes = range(1, 9)
            return ['%s%s' % (i[0], i[1]) for i in list(product(row_prefixes, col_prefixes))]
        def get_96_well_ids():
            row_prefixes = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
            col_prefixes = range(1, 13)
            return ['%s%s' % (i[0], i[1]) for i in list(product(row_prefixes, col_prefixes))]
        def get_384_well_ids():
            row_prefixes = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
            col_prefixes = range(1, 25)
            return ['%s%s' % (i[0], i[1]) for i in list(product(row_prefixes, col_prefixes))]

        curr_plate_id = 1
        if plate_size == PlateWells.WELLS_24:
            curr_well_ids = get_24_well_ids()
        elif plate_size == PlateWells.WELLS_48:
            curr_well_ids = get_48_well_ids()
        elif plate_size == PlateWells.WELLS_96:
            curr_well_ids = get_96_well_ids()
        elif plate_size == PlateWells.WELLS_384:
            curr_well_ids = get_384_well_ids()
        else:
            logger.error('Unsupported number of wells: %i' % plate_size)
            exit(1)

        # Step 3. Assign well IDs
        for pool_id in sorted(self.to_dataframe()['pool_id'].unique()):
            if len(curr_well_ids) == 0:
                if plate_size == PlateWells.WELLS_24:
                    curr_well_ids = get_24_well_ids()
                elif plate_size == PlateWells.WELLS_48:
                    curr_well_ids = get_48_well_ids()
                elif plate_size == PlateWells.WELLS_96:
                    curr_well_ids = get_96_well_ids()
                elif plate_size == PlateWells.WELLS_384:
                    curr_well_ids = get_384_well_ids()
                else:
                    logger.error('Unsupported number of wells: %i' % plate_size)
                    exit(1)
                curr_plate_id += 1
            curr_well_id = curr_well_ids[0]
            curr_well_ids.pop(0)
            self.plate_ids[pool_id] = (curr_plate_id, curr_well_id)

    def get_peptide_sequence(self, peptide_id: PeptideId) -> PoolId:
        """
        Returns the peptide sequence of a peptide ID.

        Returns
        -------
        peptide_sequence    :   Peptide sequence.
        """
        df = self.to_dataframe()
        return df.loc[df['peptide_id'] == peptide_id, 'peptide_sequence'].values[0]
    
    def get_pool_ids(self, peptide_id: PeptideId) -> PoolId:
        """
        Returns pool IDs for a peptide ID.

        Returns
        -------
        pool_ids    :   List of PooldId.
        """
        df = self.to_dataframe()
        return df.loc[df['peptide_id'] == peptide_id, 'pool_id'].values.tolist()

    def num_violations(self) -> float:
        """
        Number of violations
        (i.e. number of times two peptides appear together more than once).

        Returns
        -------
        num_violations      :   Number of violations.
        """
        # Step 1. Create a dictionary of peptides and pools
        df_assignments = self.to_dataframe()
        peptide_pool_dict = defaultdict(list)
        peptide_ids = list(df_assignments['peptide_id'].unique())
        for peptide_id in peptide_ids:
            peptide_pool_dict[peptide_id] = list(df_assignments.loc[df_assignments['peptide_id'] == peptide_id, 'pool_id'].unique())

        # Step 2. Enumerate the number of violations
        num_violations = 0
        for i in range(0, len(peptide_ids)):
            for j in range(i + 1, len(peptide_ids)):
                p1_pools = peptide_pool_dict[peptide_ids[i]]
                p2_pools = peptide_pool_dict[peptide_ids[j]]
                shared_pools = set(p1_pools).intersection(set(p2_pools))
                if len(shared_pools) > 1:
                    num_violations += len(shared_pools) - 1
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
        df = pd.DataFrame(data)
        df.sort_values(by=['peptide_id'], inplace=True)
        return df

    def to_bench_ready_dataframe(self) -> pd.DataFrame:
        data = {
            'plate_id': [],
            'well_id': [],
            'peptide_ids': [],
            'peptide_sequences': []
        }
        df_assignment = self.to_dataframe()
        for name, group in df_assignment.groupby(['plate_id', 'well_id']):
            data['plate_id'].append(name[0])
            data['well_id'].append(name[1])
            data['peptide_ids'].append(';'.join(group['peptide_id'].values.tolist()))
            data['peptide_sequences'].append(';'.join(group['peptide_sequence'].values.tolist()))
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
        3. Two peptides are not pooled together more than once.
        4. There is an optimal (minimal) number of pools.

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
                    if constraint_2_bool:
                        logger.info("Assignment does not meet constraint #2: there are peptides that do not belong to exactly one unique combination of pool IDs.")
                    logger.info("\tPools %s have the following peptides: %s." % (key, ','.join(value)))
                constraint_2_bool = False
        if constraint_2_bool:
            if verbose:
                logger.info('Assignment meets constraint #2: each peptide belongs to exactly one unique combination of pool IDs.')

        # Step 3. Two peptides are not pooled together more than once.
        constraint_3_bool = True
        num_violations = self.num_violations()
        if num_violations > 0:
            constraint_3_bool = False
            if verbose:
                logger.info("Assignment does not meet constraint #3: violation score is %i "
                            "(proxy of number of times peptide pairs are pooled together more than once)." %
                            num_violations)
        else:
            if verbose:
                logger.info('Assignment meets constraint #3: every pair of peptides is pooled together at most once.')

        # Step 4. Check that there is an optimal number of pools
        constraint_4_bool = True
        num_pools = math.ceil(len(df_assignments['peptide_id'].unique()) / num_peptides_per_pool) * num_coverage
        if len(df_assignments['pool_id'].unique()) != num_pools:
            num_extra_pools = len(df_assignments['pool_id'].unique()) - num_pools
            if verbose:
                logger.info('Assignment does not meet constraint #4: %i extra pool(s) than the minimum possible number of pools (%i).' %
                            (num_extra_pools, num_pools))
            constraint_4_bool = False
        if constraint_4_bool:
            if verbose:
                logger.info('Assignment meets constraint #4: there is an optimal (minimal) number of pools (%i).' % num_pools)

        return constraint_1_bool & constraint_2_bool & constraint_3_bool & constraint_4_bool

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

    def to_golfy_design(self) -> Tuple[Design, PeptideIndices]:
        """
        Returns in golfy design variable type.

        Returns
        -------
        design          :   Design object.
        peptide_indices :   Dictionary where
                            key     = peptide index
                            value   = peptide ID
        """
        # Step 1. Create dictionaries of peptide indices to IDs (and vice and versa)
        peptide_idx_to_id_dict = {}
        peptide_id_to_idx_dict = {}
        idx = 0
        for peptide_id in self.peptide_ids:
            peptide_idx_to_id_dict[idx] = peptide_id
            peptide_id_to_idx_dict[peptide_id] = idx
            idx += 1

        # Step 2. Create assignments for golfy Design class
        assignments = {}
        num_peptides_per_pool = -1
        for coverage in self.assignments.keys():
            assignments[coverage-1] = {}
            for pool in self.assignments[coverage].keys():
                assignments[coverage-1][pool-1] = []
                for peptide_id, peptide_sequence in self.assignments[coverage][pool]:
                    assignments[coverage-1][pool-1].append(peptide_id_to_idx_dict[peptide_id])
                if num_peptides_per_pool < len(self.assignments[coverage][pool]):
                    num_peptides_per_pool = len(self.assignments[coverage][pool])

        # Step 3. Create a golfy Design object
        design = Design(
            num_peptides=len(self.peptide_ids),
            max_peptides_per_pool=num_peptides_per_pool,
            num_replicates=len(self.assignments.keys()),
            allow_extra_pools=False,
            invalid_neighbors=[],
            preferred_neighbors=[],
            assignments=assignments
        )
        return design, peptide_idx_to_id_dict

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
        for peptide_neighbor in peptide_neighbors:
            random.shuffle(peptide_neighbor)

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
            sheet_name: str = 'assignment'
    ) -> 'BlockAssignment':
        """
        Reads an Excel file and returns a BlockAssignment object.

        Parameters
        ----------
        excel_file          :   Excel file.
        sheet_name          :   Sheet name (default: 'assignment').

        Returns
        -------
        block_assignment    :   BlockAssignment object.
        """
        df = pd.read_excel(excel_file, sheet_name=sheet_name, engine='openpyxl')

        # Step 1. Generate coverage IDs
        pool_idx = 1
        plate_well_ids = {} # key   = <plate_id>-<well_id>
                            # value = <pool_id>
        for index, row in df.iterrows():
            curr_plate_well_id = '%s-%s' % (row['plate_id'], row['well_id'])
            if curr_plate_well_id not in plate_well_ids.keys():
                plate_well_ids[curr_plate_well_id] = pool_idx
                pool_idx += 1

        # Step 2. Add peptides
        block_assignment = BlockAssignment()
        coverage_ids = {}   # key   = <peptide_id>
                            # value = count enumerated
        for index, row in df.iterrows():
            curr_peptide_id = row['peptide_id']
            curr_peptide_sequence = row['peptide_sequence']
            curr_plate_id = row['plate_id']
            curr_well_id = row['well_id']
            curr_pool_id = plate_well_ids['%s-%s' % (curr_plate_id, curr_well_id)]
            curr_coverage_id = coverage_ids.get(curr_peptide_id, 0) + 1
            coverage_ids[curr_peptide_id] = curr_coverage_id
            block_assignment.add_peptide(
                coverage=curr_coverage_id,
                pool=curr_pool_id,
                peptide_id=curr_peptide_id,
                peptide_sequence=curr_peptide_sequence
            )
            block_assignment.plate_ids[curr_pool_id] = (curr_plate_id, curr_well_id)
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
        Minimizes violations (i.e. number of times peptide pairs are pooled together more than once) in a list of
        block assignments by shuffling pool IDs.

        Parameters
        ----------
        block_assignments   :   List of BlockAssignment objects.
        shuffle_iters       :   Number of iterations to shuffle pool IDs.

        Returns
        -------
        block_assignments   :   List of BlockAssignment objects.
        """
        min_violations = BlockAssignment.merge(block_assignments=block_assignments).num_violations()
        curr_block_assignments = copy.deepcopy(block_assignments)
        best_block_assignments = copy.deepcopy(block_assignments)
        for _ in range(0, shuffle_iters):
            random_idx = random.choice(list(range(0, len(curr_block_assignments))))
            curr_block_assignments[random_idx].shuffle_pool_ids()
            curr_num_violations = BlockAssignment.merge(block_assignments=curr_block_assignments).num_violations()
            if curr_num_violations < min_violations:
                if verbose:
                    logger.info('\tFound a better assignment; current number of violations: %i, new number of violations: %i' %
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
    peptide_pairs       :   List of tuples (peptide ID, peptide ID, score).

    Returns
    -------
    peptide_neighbors   :   List of list of peptide IDs.
    """
    # Step 1. Create a dictionary to store peptide relationships
    peptide_dict = {}
    for peptide_id_1, peptide_id_2, score in peptide_pairs:
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


