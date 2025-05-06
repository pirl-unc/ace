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
import math
import random
from collections import defaultdict
from dataclasses import dataclass, field
from golfy import Design
from itertools import combinations, product
from typing import Dict, List, Mapping, Tuple
from .constants import *
from .logger import get_logger
from .peptide import Peptide
from .peptide_set import PeptideSet
from .plate_well import PlateWell
from .pool import Pool
from .types import *


logger = get_logger(__name__)


def compute_transitive_neighbors(
        peptide_pairs: List[Tuple[str,str,float]]
) -> List[List[str]]:
    """
    Computes transitive peptide neighbors.

    Parameters:
        peptide_pairs       :   List of (peptide ID, peptide ID).

    Returns:
        peptide_neighbors   :   List of list of peptide IDs.
    """
    # Step 1. Create a dictionary to store peptide relationships
    # peptide_dict:
    # {
    #     peptide_id_1: set(peptide_id_2, peptide_id_3,...),
    #     peptide_id_2: set(peptide_id_1, peptide_id_4,...),
    #     ...
    # }
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


def infer_coverage_ids(df: pd.DataFrame) -> pd.DataFrame:
    """
    Infer coverage IDs.

    Parameters:
        df      :   Pandas DataFrame with the following columns:

                        - 'peptide_id'
                        - 'peptide_sequence'
                        - 'plate_id'
                        - 'well_id'

    Returns:
        df      :   Pandas DataFrame with the following columns:

                        - 'peptide_id'
                        - 'peptide_sequence'
                        - 'plate_id'
                        - 'well_id'
                        - 'coverage_id'
    """
    df['pool'] = df['plate_id'].astype(str) + ':' + df['well_id'].astype(str)

    # Step 1: Infer the coverage
    num_coverage = 0
    for peptide_id in df['peptide_id'].unique():
        df_matched = df[df['peptide_id'] == peptide_id]
        if len(df_matched) > num_coverage:
            num_coverage = len(df_matched)
    print("Inferred coverage: %i" % num_coverage)
    for peptide_id in df['peptide_id'].unique():
        df_matched = df[df['peptide_id'] == peptide_id]
        if len(df_matched) != num_coverage:
            raise Exception('Coverage is different for %s' % peptide_id)

    # Step 2: Map peptide to pools
    peptide_pools = {}
    for peptide_id in df['peptide_id'].unique():
        pool_ids = df.loc[df['peptide_id'] == peptide_id, 'pool'].values.tolist()
        peptide_pools[peptide_id] = set(pool_ids)

    # Step 3: Map pool to peptides
    pool_peptides = {}
    for pool_id in df['pool'].unique():
        peptide_ids = df.loc[df['pool'] == pool_id, 'peptide_id'].values.tolist()
        pool_peptides[pool_id] = set(peptide_ids)

    # Step 4: Initialize pool to coverage assignments
    all_pools = list(df['pool'].unique())
    pool_to_coverage = {}
    coverage_to_pools = {i: [] for i in range(1, num_coverage + 1)}

    # Step 5: Solve using backtracking
    solution_found = solve_coverage_assignment(
        pool_idx=0,
        all_pools=all_pools,
        pool_to_coverage=pool_to_coverage,
        coverage_to_pools=coverage_to_pools,
        peptide_pools=peptide_pools,
        pool_peptides=pool_peptides,
        num_coverage=num_coverage
    )

    if solution_found:
        df['coverage_id'] = df['pool'].map(pool_to_coverage)
        return df
    else:
        raise Exception('Coverage IDs could not be assigned.')


def solve_coverage_assignment(
        pool_idx: int,
        all_pools: List[str],
        pool_to_coverage: Dict[str, int],
        coverage_to_pools: Dict[int, List[str]],
        peptide_pools: Dict[str, str],
        pool_peptides: Dict[str, str],
        num_coverage: int
):
    """
    Recursively assign coverage IDs to pools using backtracking.

    Parameters:
        pool_idx            :   Index of current pool being assigned.
        all_pools           :   List of all pool IDs.
        pool_to_coverage    :   Dictionary mapping each pool to its assigned coverage ID.
        coverage_to_pools   :   Dictionary mapping each coverage ID to its assigned pools.
        peptide_pools       :   Dictionary mapping each peptide to its pools.
        pool_peptides       :   Dictionary mapping each pool to its peptides.
        num_coverage        :   Number of coverage batches.

    Returns:
        bool: Whether a valid assignment was found
    """
    # All pools have been assigned
    if pool_idx == len(all_pools):
        # Check if each peptide appears in each coverage ID
        for peptide, pools in peptide_pools.items():
            coverage_appearances = set()
            for pool in pools:
                if pool in pool_to_coverage:
                    coverage_appearances.add(pool_to_coverage[pool])
            if len(coverage_appearances) != num_coverage:
                return False
        return True

    current_pool = all_pools[pool_idx]
    current_peptides = pool_peptides[current_pool]

    # Try each possible coverage ID for the current pool
    for coverage_id in range(1, num_coverage + 1):
        # Check if this assignment would be valid
        is_valid = is_valid_assignment(
            peptides=current_peptides,
            coverage_id=coverage_id,
            pool_to_coverage=pool_to_coverage,
            peptide_pools=peptide_pools
        )

        if is_valid:
            # Make the assignment
            pool_to_coverage[current_pool] = coverage_id
            coverage_to_pools[coverage_id].append(current_pool)

            # Recursively try to assign the next pool
            solution_found = solve_coverage_assignment(
                pool_idx=pool_idx + 1,
                all_pools=all_pools,
                pool_to_coverage=pool_to_coverage,
                coverage_to_pools=coverage_to_pools,
                peptide_pools=peptide_pools,
                pool_peptides=pool_peptides,
                num_coverage=num_coverage
            )

            if solution_found:
                return True

            # If the recursive call did not find a solution, backtrack
            pool_to_coverage.pop(current_pool)
            coverage_to_pools[coverage_id].remove(current_pool)

    # If no valid assignment was found for current_pool
    return False


def is_valid_assignment(
        peptides,
        coverage_id,
        pool_to_coverage,
        peptide_pools
) -> bool:
    """
    Check if assigning the given coverage_id to the pool is valid.

    Parameters:
        peptides            :   Set of peptides in the pool.
        coverage_id         :   The coverage ID being assigned.
        pool_to_coverage    :   Current pool to coverage assignments.
        coverage_to_pools   :   Current coverage to pools assignments.
        peptide_pools       :   Dictionary mapping each peptide to its pools.

    Returns:
        bool                :   Whether the assignment is valid.
    """
    # Check if any peptide in this pool already has a pool with this coverage_id
    for peptide in peptides:
        for p in peptide_pools[peptide]:
            if p in pool_to_coverage and pool_to_coverage[p] == coverage_id:
                return False
    return True


@dataclass
class BlockAssignment:
    pools: Mapping[int, Pool] = field(
        default_factory=dict,
        metadata={"doc": "Mapping from a pool ID to a Pool object."}
    )
    plate_map: Mapping[int,PlateWell] = field(
        default_factory=dict,
        metadata={"doc": "Mapping from a pool ID to plate and well IDs."}
    )

    @property
    def coverage_ids(self) -> List[int]:
        """
        Return all coverage IDs.

        Returns:
            coverage_ids    :   Coverage IDs.
        """
        df_assignments = self.to_dataframe()
        return list(df_assignments['coverage_id'].unique())

    @property
    def num_peptides(self) -> int:
        """
        Return the total number of peptides in this BlockAssignment.

        Returns:
            num_peptides   :   Total number of peptides.
        """
        df_assignments = self.to_dataframe()
        return len(list(df_assignments['peptide_id'].unique()))

    @property
    def num_pools(self) -> int:
        """
        Return the total number of pools in this BlockAssignment.

        Returns:
            num_pools   :   Total number of pools.
        """
        df_assignments = self.to_dataframe()
        return len(list(df_assignments['pool_id'].unique()))

    @property
    def num_violations(self) -> int:
        """
        Return the number of violations (defined as number of times two peptides appear together more than once).

        Returns:
            num_violations  :   Number of violations.
        """
        # Step 1. Create a dictionary of peptides and pools
        df_assignments = self.to_dataframe()
        peptide_pool_dict = defaultdict(list) # key = peptide ID, value = list of pool IDs
        peptide_ids = list(df_assignments['peptide_id'].unique())
        for peptide_id in peptide_ids:
            peptide_pool_dict[peptide_id] = list(df_assignments.loc[df_assignments['peptide_id'] == peptide_id, 'pool_id'].unique())

        # Step 2. Count the number of violations
        num_violations = 0
        for i in range(0, len(peptide_ids)):
            for j in range(i + 1, len(peptide_ids)):
                p1_pools = peptide_pool_dict[peptide_ids[i]]
                p2_pools = peptide_pool_dict[peptide_ids[j]]
                shared_pools = set(p1_pools).intersection(set(p2_pools))
                if len(shared_pools) > 1:
                    num_violations += len(shared_pools) - 1

        return num_violations

    @property
    def peptide_ids(self) -> List[str]:
        """
        Return all peptide IDs.

        Returns:
            peptide_ids :   Peptide IDs.
        """
        df_assignments = self.to_dataframe()
        return list(df_assignments['peptide_id'].unique())

    @property
    def pool_ids(self) -> List[int]:
        """
        Return all pool IDs.

        Returns:
            pool_ids    :   Pool IDs.
        """
        df_assignments = self.to_dataframe()
        return list(df_assignments['pool_id'].unique())

    @property
    def pooled_peptide_pairs(self) -> List[Tuple[str,str]]:
        """
        Return a list of peptide pairs that are pooled together.

        Returns:
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
            peptide_id: str,
            peptide_sequence: str,
            coverage_id: int,
            pool_id: int
    ):
        """
        Add a peptide.

        Parameters:
            peptide_id          :   Peptide ID.
            peptide_sequence    :   Peptide sequence.
            coverage_id         :   Coverage ID.
            pool_id             :   Pool ID.
        """
        peptide = Peptide(id=peptide_id, sequence=peptide_sequence)
        if pool_id not in self.pools.keys():
            pool = Pool(
                id=pool_id,
                coverage_id=coverage_id
            )
            pool.add_peptide(peptide=peptide)
            self.pools[pool_id] = pool
        else:
            assert self.pools[pool_id].coverage_id == coverage_id, (
                "Coverage ID mismatch for pool ID %i: expected %i, found %i" % (pool_id, self.pools[pool_id].coverage_id, coverage_id)
            )
            self.pools[pool_id].add_peptide(peptide=peptide)

    def assign_well_ids(self, num_plate_wells: NumPlateWells):
        """
        Assign plate and well IDs to an ELISpot configuration.

        Parameters:
            num_plate_wells     :   Number of plate wells (allowed values: 24, 48, 96, 384).
        """
        # Step 1. Clear the current self.plate_ids
        self.plate_map = {} # key = pool ID, value = (plate ID, well ID)

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
        if num_plate_wells == NumPlateWells.WELLS_24:
            curr_well_ids = get_24_well_ids()
        elif num_plate_wells == NumPlateWells.WELLS_48:
            curr_well_ids = get_48_well_ids()
        elif num_plate_wells == NumPlateWells.WELLS_96:
            curr_well_ids = get_96_well_ids()
        elif num_plate_wells == NumPlateWells.WELLS_384:
            curr_well_ids = get_384_well_ids()
        else:
            raise Exception("Unsupported number of wells: %i" % num_plate_wells)

        # Step 3. Assign well IDs
        for pool_id in sorted(self.to_dataframe()['pool_id'].unique()):
            if len(curr_well_ids) == 0:
                if num_plate_wells == NumPlateWells.WELLS_24:
                    curr_well_ids = get_24_well_ids()
                elif num_plate_wells == NumPlateWells.WELLS_48:
                    curr_well_ids = get_48_well_ids()
                elif num_plate_wells == NumPlateWells.WELLS_96:
                    curr_well_ids = get_96_well_ids()
                elif num_plate_wells == NumPlateWells.WELLS_384:
                    curr_well_ids = get_384_well_ids()
                else:
                    raise Exception("Unsupported number of wells: %i" % num_plate_wells)
                curr_plate_id += 1
            curr_well_id = curr_well_ids[0]
            curr_well_ids.pop(0)
            plate_well = PlateWell(
                plate_id=curr_plate_id,
                well_id=curr_well_id
            )
            self.plate_map[pool_id] = plate_well

    def get_peptide_sequence(self, peptide_id: str) -> str:
        """
        Return the peptide sequence given a peptide ID.

        Returns:
            peptide_sequence    :   Peptide sequence.
        """
        df = self.to_dataframe()
        return str(df.loc[df['peptide_id'] == peptide_id, 'peptide_sequence'].values[0])

    def get_pool_ids(self, peptide_id: str) -> List[int]:
        """
        Return pool IDs for a peptide ID.

        Returns:
            pool_ids    :   List of pool IDs.
        """
        df = self.to_dataframe()
        return df.loc[df['peptide_id'] == peptide_id, 'pool_id'].values.tolist()

    def is_optimal(
            self,
            num_coverage: int,
            num_peptides_per_pool: int,
            verbose: bool = True
    ) -> bool:
        """
        Verify whether a given ELISpot assignment satisfies the following constraints:
            1. Each peptide is in 'num_coverage' number of different pools.
            2. Each peptide is in exactly one unique combination of pool IDs.
            3. Two peptides are not pooled together more than once.
            4. There is an optimal (minimal) number of pools.

        Parameters:
            num_coverage            :   Coverage.
            num_peptides_per_pool   :   Number of peptides per pool.
            verbose                 :   Verbose.

        Returns:
            is_optimal              :   True if the input configuration meets all desired criteria. False otherwise.
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
        num_violations = self.num_violations
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
        Shuffle pool IDs.
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
        old_pools = self.pools
        self.pools = {}
        for old_pool_id, pool in old_pools.items():
            new_pool_id = new_pool_ids_dict[old_pool_id]
            coverage_id = pool.coverage_id
            for peptide in pool.peptides:
                self.add_peptide(
                    peptide_id=peptide.id,
                    peptide_sequence=peptide.sequence,
                    coverage_id=coverage_id,
                    pool_id=new_pool_id
                )

    def load_plate_map(self, plate_map: Dict[int,PlateWell]):
        self.plate_map = plate_map

    def to_bench_ready_dataframe(self) -> pd.DataFrame:
        """
        Return a bench-ready Pandas DataFrame of this BlockAssignment.

        Returns:
            df_assignments  :   Pandas DataFrame with the following columns:

                                    - 'plate_id'
                                    - 'well_id'
                                    - 'peptide_ids'
                                    - 'peptide_sequences'
        """
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

    def to_dataframe(self) -> pd.DataFrame:
        """
        Return a Pandas DataFrame of this BlockAssignment.

        Returns:
            df_assignments  :   Pandas DataFrame with the following columns:

                                    - 'coverage_id'
                                    - 'pool_id'
                                    - 'peptide_id'
                                    - 'peptide_sequence'
                                    - 'plate_id'
                                    - 'well_id'
        """
        data = {
            'coverage_id': [],
            'pool_id': [],
            'peptide_id': [],
            'peptide_sequence': [],
            'plate_id': [],
            'well_id': []
        }
        for pool_id,pool in self.pools.items():
            coverage_id = pool.coverage_id
            pool_id = pool.id
            for peptide in pool.peptides:
                data['coverage_id'].append(coverage_id)
                data['pool_id'].append(pool_id)
                data['peptide_id'].append(peptide.id)
                data['peptide_sequence'].append(peptide.sequence)
                if pool_id in self.plate_map.keys():
                    data['plate_id'].append(self.plate_map[pool_id].plate_id)
                    data['well_id'].append(self.plate_map[pool_id].well_id)
                else:
                    data['plate_id'].append('')
                    data['well_id'].append('')
        df = pd.DataFrame(data)
        df.sort_values(by=['peptide_id'], inplace=True)
        return df

    def to_golfy_design(self) -> Tuple[Design, Mapping[int,str]]:
        """
        Convert to a Golfy design representation.

        Returns:
            Tuple[Design, PeptideIndices]:
                - A Design object compatible with the golfy package.
                - A dictionary mapping peptide indices to peptide IDs.
        """
        # Step 1. Create dictionaries of peptide indices to IDs (and vice and versa)
        peptide_idx_to_id_dict = {}
        peptide_id_to_idx_dict = {}
        for idx, peptide_id in enumerate(self.peptide_ids):
            peptide_idx_to_id_dict[idx] = peptide_id
            peptide_id_to_idx_dict[peptide_id] = idx

        # Step 2. Create assignments for golfy Design class
        # assignments:
        # {
        #     coverage_id: {
        #         pool_id: [
        #             peptide_index_1,
        #             peptide_index_2,
        #             ...
        #         ]
        #     }
        # }
        assignments = {}
        coverage_ids = set()
        num_peptides_per_pool = -1
        for pool_id, pool in self.pools.items():
            coverage_id = pool.coverage_id
            coverage_ids.add(coverage_id)

            # Zero-indexed coverage and pool
            cidx = coverage_id - 1
            pidx = pool_id - 1

            if cidx not in assignments:
                assignments[cidx] = {}
            assignments[cidx][pidx] = [peptide_id_to_idx_dict[pep.id] for pep in pool.peptides]

            num_peptides_per_pool = max(num_peptides_per_pool, pool.num_peptides)

        # Step 3. Create a golfy Design object
        design = Design(
            num_peptides=self.num_peptides,
            max_peptides_per_pool=num_peptides_per_pool,
            num_replicates=len(coverage_ids),
            allow_extra_pools=False,
            invalid_neighbors=[],
            preferred_neighbors=[],
            assignments=assignments
        )
        return design, peptide_idx_to_id_dict

    @staticmethod
    def generate_single_coverage_block_assignment(
            peptides: List[Peptide],
            preferred_peptide_pairs: List[Tuple[str,str,float]],
            num_peptides_per_pool: int,
            coverage: int = 1
    ) -> 'BlockAssignment':
        """
        Generate a one-coverage block assignment given a list of preferred peptide pairs.

        Parameters:
            peptides                    :   List of Peptide objects.
            preferred_peptide_pairs     :   List of preferred peptide pairs (peptide IDs).
            num_peptides_per_pool       :   Number of peptides per pool.
            coverage                    :   Coverage (default: 1).

        Returns:
            block_assignment            :   BlockAssignment object.
        """
        peptide_set = PeptideSet()
        for peptide in peptides:
            peptide_set.add_peptide(peptide=peptide)
        df_peptides = peptide_set.to_dataframe()

        # Step 1. Initialize pools
        # pools:
        # {
        #   pool_id: [
        #       peptide_id_1,
        #       peptide_id_2,
        #       ...
        #   ]
        # }
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
        for pool_id, peptide_ids in pools.items():
            for peptide_id in peptide_ids:
                peptide_sequence = df_peptides.loc[df_peptides['peptide_id'] == peptide_id, 'peptide_sequence'].values[
                    0]
                block_assignment.add_peptide(
                    peptide_id=peptide_id,
                    peptide_sequence=peptide_sequence,
                    coverage_id=coverage,
                    pool_id=pool_id
                )
        return block_assignment

    @staticmethod
    def load_from_dataframe(df_assignment: pd.DataFrame) -> 'BlockAssignment':
        """
        Load block assignments from a Pandas DataFrame and return a BlockAssignment object.

        Parameters:
            df_assignment       :   Pandas DataFrame with the following columns:

                                        - 'coverage_id'
                                        - 'pool_id'
                                        - 'peptide_id'
                                        - 'peptide_sequence'

        Returns:
            block_assignment    :   BlockAssignment object.
        """
        block_assignment = BlockAssignment()
        for index, row in df_assignment.iterrows():
            coverage_id = int(row['coverage_id'])
            pool_id = int(row['pool_id'])
            peptide_id = str(row['peptide_id'])
            peptide_sequence = str(row['peptide_sequence'])
            block_assignment.add_peptide(
                peptide_id=peptide_id,
                peptide_sequence=peptide_sequence,
                coverage_id=coverage_id,
                pool_id=pool_id
            )
            if 'plate_id' in row and 'well_id' in row:
                plate_id = row['plate_id']
                well_id = row['well_id']
                block_assignment.plate_map[pool_id] = PlateWell(plate_id=plate_id, well_id=well_id)
        return block_assignment

    @staticmethod
    def load_golfy_assignment(
            golfy_assignment: Dict,
            df_peptides: pd.DataFrame
    ) -> 'BlockAssignment':
        """
        Convert golfy results into a Pandas DataFrame.

        Parameters:
            golfy_assignment    :   Dictionary.
            df_peptides         :   Pandas DatFrame with the following columns:

                                        - 'peptide_id'
                                        - 'peptide_index'
                                        - 'peptide_sequence'

        Returns:
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
                        peptide_id=curr_peptide_id,
                        peptide_sequence=curr_peptide_sequence,
                        coverage_id=curr_coverage_id,
                        pool_id=curr_pool_id
                    )
                curr_pool_id += 1
        return block_assignment

    @staticmethod
    def read_excel_file(
            excel_file: str,
            sheet_name: str = 'assignment'
    ) -> 'BlockAssignment':
        """
        Read an Excel file and return a BlockAssignment object.

        Parameters:
            excel_file          :   Excel file.
            sheet_name          :   Sheet name (default: 'assignment').

        Returns:
            block_assignment    :   BlockAssignment object.
        """
        df = pd.read_excel(excel_file, sheet_name=sheet_name, engine='openpyxl')

        # Step 1. Get coverage IDs
        if 'coverage_id' not in df.columns:
            df = infer_coverage_ids(df=df)
        plate_well_coverages = {}
        for _,row in df.iterrows():
            curr_plate_id = row['plate_id']
            curr_well_id = row['well_id']
            curr_coverage_id = row['coverage_id']
            plate_well_key = '%i-%s' % (curr_plate_id, curr_well_id)
            if plate_well_key not in plate_well_coverages:
                plate_well_coverages[plate_well_key] = curr_coverage_id

        # Step 1. Assign pool IDs
        # plate_well_pools:
        # {
        #     <plate_id>-<well_id>: pool_id_1,
        #     ...
        # }
        plate_well_pools = {}
        pool_id = 1
        for index, row in df.iterrows():
            plate_well_key = '%i-%s' % (row['plate_id'], row['well_id'])
            if plate_well_key not in plate_well_pools:
                plate_well_pools[plate_well_key] = pool_id
                pool_id += 1

        # Step 2. Add peptides
        block_assignment = BlockAssignment()
        # coverage_ids:
        # {
        #     peptide_id_1: 1,
        #     peptide_id_2: 3,
        #     ...
        # }
        for index, row in df.iterrows():
            curr_peptide_id = row['peptide_id']
            curr_peptide_sequence = row['peptide_sequence']
            curr_plate_id = row['plate_id']
            curr_well_id = row['well_id']
            plate_well_key = '%i-%s' % (row['plate_id'], row['well_id'])
            curr_pool_id = int(plate_well_pools[plate_well_key])
            curr_coverage_id = int(plate_well_coverages[plate_well_key])
            block_assignment.add_peptide(
                peptide_id=curr_peptide_id,
                peptide_sequence=curr_peptide_sequence,
                coverage_id=curr_coverage_id,
                pool_id=curr_pool_id
            )
            block_assignment.plate_map[curr_pool_id] = PlateWell(plate_id=curr_plate_id, well_id=curr_well_id)
        return block_assignment

    @staticmethod
    def update_ids(
            block_assignment: 'BlockAssignment',
            start_pool_num: int,
            start_coverage_num: int
    ) -> 'BlockAssignment':
        """
        Update pool and coverage IDs of a BlockAssignment object.

        Parameters:
            block_assignment        :   BlockAssignment object.
            start_pool_num          :   Start pool number (e.g. 1).
            start_coverage_num      :   Start coverage number (e.g. 1).

        Returns:
            new_block_assignment    :   BlockAssignment object.
        """
        # Step 1. Get old coverage and pool IDs
        old_coverage_ids = sorted(block_assignment.coverage_ids)
        old_pool_ids = sorted(block_assignment.pool_ids)

        # Step 2. Create a dictionary of old and new coverage IDs
        # coverage_ids_dict:
        # {
        #     old_coverage_id_1: new_coverage_id_1,
        #     old_coverage_id_2: new_coverage_id_2,
        #     ...
        # }
        coverage_ids_dict = {}
        curr_coverage_num = start_coverage_num
        for coverage_id in old_coverage_ids:
            coverage_ids_dict[coverage_id] = curr_coverage_num
            curr_coverage_num += 1

        # Step 3. Create a dictionary of old and new pool IDs
        # pool_ids_dict:
        # {
        #     old_pool_id_1: new_pool_id_1,
        #     old_pool_id_2: new_pool_id_2,
        #     ...
        # }
        pool_ids_dict = {}
        curr_pool_num = start_pool_num
        for pool_id in old_pool_ids:
            pool_ids_dict[pool_id] = curr_pool_num
            curr_pool_num += 1

        # Step 4. Create new assignments
        new_block_assignment = BlockAssignment()
        for pool_id,pool in block_assignment.pools.items():
            new_coverage_id = coverage_ids_dict[pool.coverage_id]
            new_pool_id = pool_ids_dict[pool.id]
            for peptide in pool.peptides:
                new_block_assignment.add_peptide(
                    peptide_id=peptide.id,
                    peptide_sequence=peptide.sequence,
                    coverage_id=new_coverage_id,
                    pool_id=new_pool_id
                )

        return new_block_assignment

    @staticmethod
    def merge(block_assignments: List['BlockAssignment']) -> 'BlockAssignment':
        """
        Merges a list of BlockAssignment objects.

        Parameters:
            block_assignments   :   List of BlockAssignment objects.

        Returns:
            block_assignment    :   BlockAssignment object.
        """
        merged_block_assignment = BlockAssignment()
        for block_assignment in block_assignments:
            for pool_id,pool in block_assignment.pools.items():
                for peptide in pool.peptides:
                    merged_block_assignment.add_peptide(
                        peptide_id=peptide.id,
                        peptide_sequence=peptide.sequence,
                        coverage_id=pool.coverage_id,
                        pool_id=pool_id
                    )
        return merged_block_assignment

    @staticmethod
    def minimize_violations(
            block_assignments: List['BlockAssignment'],
            shuffle_iters: int,
            verbose: bool = True
    ) -> List['BlockAssignment']:
        """
        Minimize violations (i.e. number of times peptide pairs are pooled together more than once) in a list of
        block assignments by shuffling pool IDs.

        Parameters:
            block_assignments   :   List of BlockAssignment objects.
            shuffle_iters       :   Number of iterations to shuffle pool IDs.

        Returns:
            block_assignments   :   List of BlockAssignment objects.
        """
        min_violations = BlockAssignment.merge(block_assignments=block_assignments).num_violations
        curr_block_assignments = copy.deepcopy(block_assignments)
        best_block_assignments = copy.deepcopy(block_assignments)
        for _ in range(0, shuffle_iters):
            random_idx = random.choice(list(range(0, len(curr_block_assignments))))
            curr_block_assignments[random_idx].shuffle_pool_ids()
            curr_num_violations = BlockAssignment.merge(block_assignments=curr_block_assignments).num_violations
            if curr_num_violations < min_violations:
                if verbose:
                    logger.info('\tFound a better assignment; current number of violations: %i, new number of violations: %i' %
                                (min_violations, curr_num_violations))
                best_block_assignments = copy.deepcopy(curr_block_assignments)
                min_violations = curr_num_violations
        return best_block_assignments
