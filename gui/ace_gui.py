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


import eel
import multiprocessing
import pandas as pd
import random
import socket
import os
import torch
from importlib import resources
from typing import Dict
from acelib.main import run_ace_generate, run_ace_deconvolve
from acelib.block_assignment import BlockAssignment
from acelib.utilities import get_open_port


eel.init('views')


@eel.expose
def generate_configuration(sequences_available: bool,
                           peptides: Dict,
                           num_peptides: int,
                           num_peptides_per_pool: int,
                           num_coverage: int,
                           plate_size: str,
                           cluster_peptides: str,
                           sequence_similarity_fxn: str,
                           sequence_similarity_threshold: float,
                           max_iters: int,
                           init_strategy: str,
                           allow_extra_pools: str):
    # Step 1. Preprocess input parameters
    num_peptides = int(num_peptides)
    peptides_ = []
    if sequences_available:
        for peptide in peptides:
            peptides_.append((str(peptide['id']), peptide['sequence']))
    else:
        peptides_ = []
        for i in range(0, num_peptides):
            peptides_.append((str(i+1), ''))
    num_peptides_per_pool = int(num_peptides_per_pool)
    num_coverage = int(num_coverage)
    if sequences_available and cluster_peptides == 'yes':
        cluster_peptides = True
    else:
        cluster_peptides = False
    sequence_similarity_threshold = float(sequence_similarity_threshold)
    max_iters = int(max_iters)
    if allow_extra_pools == 'yes':
        allow_extra_pools = True
    else:
        allow_extra_pools = False

    # Step 2. Generate an ELISpot configuration
    # trained_model_file = resources.path('acelib.resources.models', 'trained_model_w_data_augmentation_b3000.pt'),
    dir_path = os.path.dirname(os.path.realpath(__file__))
    print(dir_path)
    block_assignment, block_design = run_ace_generate(
        peptides=peptides_,
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage,
        cluster_peptides=cluster_peptides,
        trained_model_file=dir_path + '/trained_model_w_data_augmentation_b3000.pt',
        golfy_allow_extra_pools=allow_extra_pools,
        sequence_similarity_function=sequence_similarity_fxn,
        sequence_similarity_threshold=sequence_similarity_threshold,
        mode='golfy',
        golfy_random_seed=random.randint(1,100000000),
        golfy_max_iters=max_iters,
        golfy_strategy=init_strategy,
        plate_size=plate_size
    )

    return block_assignment.to_dataframe().to_dict(), \
           block_assignment.to_bench_ready_dataframe().to_dict(), \
           block_design.peptides_dataframe.to_dict(), \
           block_design.metadata_dataframe.to_dict(),\
           block_design.preferred_peptide_pairs_dataframe.to_dict()


@eel.expose
def deconvolve(assignments,
               readouts,
               statistical_deconvolution_method,
               min_coverage,
               min_spot_count):
    """
    Deconvolves positive peptides.

    Parameters
    ----------
    assignments                         :   Dictionary of lists with the following keys:
                                            'peptide_id'
                                            'peptide_sequence'
                                            'plate_id'
                                            'well_id'
    readouts                            :   Dictionary of lists with the following keys:
                                            'plate_id'
                                            'well_id'
                                            'spot_count'
    statistical_deconvolution_method
    min_coverage
    min_spot_count

    Returns
    -------

    """
    # Step 1. Process assignments
    data = {
        'peptide_id': [],
        'peptide_sequence': [],
        'plate_id': [],
        'well_id': []
    }
    for assignment in assignments:
        data['peptide_id'].append(assignment['peptide_id'])
        data['peptide_sequence'].append(assignment['peptide_sequence'])
        data['plate_id'].append(assignment['plate_id'])
        data['well_id'].append(assignment['well_id'])
    df_assignment = pd.DataFrame(data)
    pool_ids_dict = {}  # key   = <plate_id>-<well_id>
                        # value = <pool_id>
    well_ids_dict = {}  # key   = <pool_id>
                        # value = <plate_id>-<well_id>
    pool_idx = 1
    for name, _ in df_assignment.groupby(['plate_id', 'well_id']):
        curr_key = '%s-%s' % (name[0], name[1])
        if curr_key not in pool_ids_dict.keys():
            pool_ids_dict[curr_key] = pool_idx
            well_ids_dict[pool_idx] = curr_key
            pool_idx += 1
    peptide_counts = {} # key   = <peptide_id>
                        # value = count
    pool_ids = []
    coverage_ids = []
    for index, row in df_assignment.iterrows():
        coverage_id = peptide_counts.get(row['peptide_id'], 1)
        peptide_counts[row['peptide_id']] = coverage_id + 1
        pool_id = pool_ids_dict['%s-%s' % (row['plate_id'], row['well_id'])]
        pool_ids.append(pool_id)
        coverage_ids.append(coverage_id)
    df_assignment['pool_id'] = pool_ids
    df_assignment['coverage_id'] = coverage_ids
    block_assignment = BlockAssignment()
    for index, row in df_assignment.iterrows():
        block_assignment.add_peptide(
            coverage=row['coverage_id'],
            pool=row['pool_id'],
            peptide_id=row['peptide_id'],
            peptide_sequence=row['peptide_sequence']
        )

    # Step 2. Process readout
    data = {
        'plate_id': [],
        'well_id': [],
        'pool_id': [],
        'spot_count': []
    }
    for readout in readouts:
        data['plate_id'].append(readout['plate_id'])
        data['well_id'].append(readout['well_id'])
        data['spot_count'].append(float(readout['spot_count']))
        data['pool_id'].append(pool_ids_dict['%s-%s' % (readout['plate_id'], readout['well_id'])])
    df_readout = pd.DataFrame(data)

    # Step 3. Perform empirical deconvolution
    min_coverage = int(min_coverage)
    min_spot_count = float(min_spot_count)
    deconvolution_result = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        statistical_deconvolution_method=statistical_deconvolution_method,
        statistical_min_peptide_activity=1.0,
        empirical_min_coverage=min_coverage,
        empirical_min_spot_count=min_spot_count
    )

    # Step 4. Convert pool IDs to well IDs
    data = {
        'peptide_id': [],
        'peptide_sequence': [],
        'estimated_peptide_spot_count': [],
        'hit_well_ids': [],
        'hit_well_ids_count': [],
        'deconvolution_result': []
    }
    for index, row in deconvolution_result.to_dataframe().iterrows():
        well_ids = []
        for curr_pool_id in row['hit_pool_ids'].split(';'):
            if curr_pool_id != '':
                well_ids.append(well_ids_dict[int(curr_pool_id)])
        data['peptide_id'].append(row['peptide_id'])
        data['peptide_sequence'].append(row['peptide_sequence'])
        data['estimated_peptide_spot_count'].append(row['estimated_peptide_spot_count'])
        data['hit_well_ids'].append(';'.join(well_ids))
        data['hit_well_ids_count'].append(len(well_ids))
        data['deconvolution_result'].append(row['deconvolution_result'])
    df_results = pd.DataFrame(data)
    return df_results.to_dict()


def find_port(port=1111):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        if s.connect_ex(("localhost", port)) == 0:
            return find_port(port=port + 1)
        else:
            return port


if __name__ == '__main__':
    multiprocessing.freeze_support()
    eel.start('index.html', size=(1920, 1080), port=get_open_port())
