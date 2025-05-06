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
from acelib.constants import *
from acelib.main import run_ace_generate, run_ace_deconvolve
from acelib.block_assignment import BlockAssignment, infer_coverage_ids
from acelib.plate_well import PlateWell
from acelib.peptide import Peptide
from acelib.utilities import get_open_port


eel.init('views')


@eel.expose
def generate_configuration(
        sequences_available: bool,
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
        allow_extra_pools: str
):
    # Step 1. Preprocess input parameters
    num_peptides = int(num_peptides)
    peptides_ = []
    if sequences_available:
        for peptide in peptides:
            peptide_ = Peptide(id=str(peptide['id']), sequence=str(peptide['sequence']))
            peptides_.append(peptide_)
    else:
        for i in range(0, num_peptides):
            peptide_ = Peptide(id=str(i+1), sequence='')
            peptides_.append(peptide_)
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
    dir_path = os.path.dirname(os.path.realpath(__file__))
    block_assignment, block_design = run_ace_generate(
        peptides=peptides_,
        num_peptides_per_pool=num_peptides_per_pool,
        num_coverage=num_coverage,
        trained_model_file=dir_path + '/trained_model_w_data_augmentation_b3000.pt',
        cluster_peptides=cluster_peptides,
        mode=GenerateMode.GOLFY,
        sequence_similarity_function=SequenceSimilarityFunction(sequence_similarity_fxn),
        sequence_similarity_threshold=sequence_similarity_threshold,
        golfy_random_seed=random.randint(1, 100000000),
        golfy_strategy=GolfyStrategy(init_strategy),
        golfy_max_iters=max_iters,
        golfy_allow_extra_pools=allow_extra_pools,
        num_plate_wells=NumPlateWells(plate_size)
    )

    return block_assignment.to_dataframe().to_dict(), \
           block_assignment.to_bench_ready_dataframe().to_dict(), \
           block_design.peptides_dataframe.to_dict(), \
           block_design.metadata_dataframe.to_dict(),\
           block_design.preferred_peptide_pairs_dataframe.to_dict()


@eel.expose
def deconvolve(
        assignments,
        readouts,
        statistical_deconvolution_method,
        min_coverage,
        min_spot_count
):
    for assignment in assignments:
        assert 'peptide_id' in assignment.keys()
        assert 'peptide_sequence' in assignment.keys()
        assert 'plate_id' in assignment.keys()
        assert 'well_id' in assignment.keys()

    for readout in readouts:
        assert 'plate_id' in readout.keys()
        assert 'well_id' in readout.keys()
        assert 'spot_count' in readout.keys()

    # Step 1. Process assignments
    data = {
        'peptide_id': [],
        'peptide_sequence': [],
        'plate_id': [],
        'well_id': []
    }
    for assignment in assignments:
        data['peptide_id'].append(str(assignment['peptide_id']))
        data['peptide_sequence'].append(str(assignment['peptide_sequence']))
        data['plate_id'].append(int(assignment['plate_id']))
        data['well_id'].append(str(assignment['well_id']) )
    df = pd.DataFrame(data)

    # Step 2. Assign coverage IDs
    df = infer_coverage_ids(df=df)

    # Step 3. Assign pool IDs
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

    # Step 4. Create a BlockAssignment object
    block_assignment = BlockAssignment()
    for index, row in df.iterrows():
        curr_plate_id = row['plate_id']
        curr_well_id = row['well_id']
        plate_well_key = '%i-%s' % (curr_plate_id, curr_well_id)
        pool_id = plate_well_pools[plate_well_key]
        coverage_id = row['coverage_id']
        block_assignment.add_peptide(
            peptide_id=row['peptide_id'],
            peptide_sequence=row['peptide_sequence'],
            pool_id=pool_id,
            coverage_id=coverage_id
        )
    plate_map = {}
    for key,pool_id in plate_well_pools.items():
        plate_id = int(key.split('-')[0])
        well_id = str(key.split('-')[1])
        plate_well = PlateWell(plate_id=plate_id, well_id=well_id)
        plate_map[pool_id] = plate_well
    block_assignment.load_plate_map(plate_map=plate_map)

    # Step 5. Process readout
    data = {
        'plate_id': [],
        'well_id': [],
        'pool_id': [],
        'spot_count': []
    }
    for readout in readouts:
        curr_plate_id = int(readout['plate_id'])
        curr_well_id = str(readout['well_id'])
        plate_well_key = '%i-%s' % (curr_plate_id, curr_well_id)
        pool_id = plate_well_pools[plate_well_key]
        data['plate_id'].append(readout['plate_id'])
        data['well_id'].append(readout['well_id'])
        data['spot_count'].append(float(readout['spot_count']))
        data['pool_id'].append(pool_id)
    df_readout = pd.DataFrame(data)
    print(df_readout.head(n=5))

    # Step 6. Perform empirical deconvolution
    min_coverage = int(min_coverage)
    min_spot_count = float(min_spot_count)
    deconvolved_peptide_set = run_ace_deconvolve(
        df_readout=df_readout,
        block_assignment=block_assignment,
        method=DeconvolutionMethod(statistical_deconvolution_method),
        min_coverage=min_coverage,
        min_pool_spot_count=min_spot_count
    )
    return deconvolved_peptide_set.to_dataframe().to_dict()


def find_port(port=1111):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        if s.connect_ex(("localhost", port)) == 0:
            return find_port(port=port + 1)
        else:
            return port


if __name__ == '__main__':
    multiprocessing.freeze_support()
    os.environ["KMP_DUPLICATE_LIB_OK"] = "True" # For Windows
    eel.start('html/index.html', size=(1920, 1080), port=get_open_port())
