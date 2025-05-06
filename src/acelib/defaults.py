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
DEFAULT_GENERATE_MODE = 'golfy'
DEFAULT_GENERATE_NUM_PLATE_WELLS = 96
DEFAULT_GENERATE_CLUSTER_PEPTIDES = True
DEFAULT_GENERATE_SEQUENCE_SIMILARITY_FUNCTION = 'euclidean'
DEFAULT_GENERATE_SEQUENCE_SIMILARITY_THRESHOLD = 0.7
DEFAULT_GENERATE_GOLFY_RANDOM_SEED = 42
DEFAULT_GENERATE_GOLFY_MAX_ITERS = 2000
DEFAULT_GENERATE_GOLFY_STRATEGY = 'greedy'
DEFAULT_GENERATE_GOLFY_ALLOW_EXTRA_POOLS = False
DEFAULT_GENERATE_CPSAT_SOLVER_NUM_PROCESSES = 2
DEFAULT_GENERATE_CPSAT_SOLVER_SHUFFLE_ITERS = 1000
DEFAULT_GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_BLOCK = 100
DEFAULT_GENERATE_CPSAT_SOLVER_MAX_PEPTIDES_PER_POOL = 10


"""deconvolve"""
DEFAULT_DECONVOLVE_METHOD = 'cem'
