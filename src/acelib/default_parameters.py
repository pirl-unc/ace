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
The purpose of this python3 script is to define ACE default parameters.
"""


"""generate"""
GENERATE_GOLFY_RANDOM_SEED = 42
GENERATE_GOLFY_MAX_ITERS = 2000
GENERATE_GOLFY_STRATEGY = 'greedy'
GENERATE_GOLFY_ALLOW_EXTRA_POOLS = False
GENERATE_SEQUENCE_SIMILARITY_THRESHOLD = 0.7
GENERATE_SEQUENCE_SIMILARITY_FUNCTION = 'euclidean'
GENERATE_CPSAT_SOLVER_NUM_PROCESSES = 2
GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_BLOCK = 100
GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_POOL = 10
GENERATE_CPSAT_SOLVER_SHUFFLE_ITERS = 1000
GENERATE_ASSIGN_WELL_IDS = True
GENERATE_NUM_PLATE_WELLS = 96


"""deconvolve"""
DECONVOLVE_METHOD = 'cem'
DECONVOLVE_MIN_PEPTIDE_SPOT_COUNT = 1.0
