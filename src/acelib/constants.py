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
The purpose of this python3 script is to define ACE constants.
"""


import pandas as pd
from enum import Enum, IntEnum


_AMINO_ACID_PROPERTIES = pd.DataFrame(
    map(lambda x: x.split(","),
        ("1.29,0.9,0,0.049,1.8,0,0,0.047,0.065,0.78,67,1,0,0,1;"
         "1.11,0.74,0,0.02,2.5,-2,0,0.015,0.015,0.8,86,1,1,-1,0;"
         "1.04,0.72,-1,0.051,-3.5,-2,1,0.071,0.074,1.41,91,1,0,1;"
         "1.44,0.75,-1,0.051,-3.5,-2,1,0.094,0.089,1,109,1,0,1,0;"
         "1.07,1.32,0,0.051,2.8,0,0,0.021,0.029,0.58,135,1,1,-1,0;"
         "0.56,0.92,0,0.06,-0.4,0,0,0.071,0.07,1.64,48,1,0,1,1;"
         "1.22,1.08,0,0.034,-3.2,1,1,0.022,0.025,0.69,118,1,0,-1,0;"
         "0.97,1.45,0,0.047,4.5,0,0,0.032,0.035,0.51,124,1,1,-1,0;"
         "1.23,0.77,1,0.05,-3.9,2,1,0.105,0.08,0.96,135,1,0,1,0;"
         "1.3,1.02,0,0.078,3.8,0,0,0.052,0.063,0.59,124,1,1,-1,1;"
         "1.47,0.97,0,0.027,1.9,0,0,0.017,0.016,0.39,124,1,1,1,0;"
         "0.9,0.76,0,0.058,-3.5,0,1,0.062,0.053,1.28,96,1,0,1,1;"
         "0.52,0.64,0,0.051,-1.6,0,0,0.052,0.054,1.91,90,1,0,1,0;"
         "1.27,0.8,0,0.051,-3.5,1,1,0.053,0.051,0.97,114,1,0,1,0;"
         "0.96,0.99,1,0.066,-4.5,2,1,0.068,0.059,0.88,148,1,0,1,1;"
         "0.82,0.95,0,0.057,-0.8,-1,1,0.072,0.071,1.33,73,1,0,1,1;"
         "0.82,1.21,0,0.064,-0.7,-1,0,0.064,0.065,1.03,93,1,0,0,1;"
         "0.91,1.49,0,0.049,4.2,0,0,0.048,0.048,0.47,105,1,1,-1,0;"
         "0.99,1.14,0,0.022,-0.9,1,1,0.007,0.012,0.75,163,1,1,-1,0;"
         "0.72,1.25,0,0.07,-1.3,-1,1,0.032,0.033,1.05,141,1,1,-1,1").split(";")),
    columns=['alpha', 'beta', 'charge', 'core', 'hydropathy', 'pH', 'polarity', 'rim', 'surface', 'turn', 'volume', 'strength', 'disorder', 'high_contact', 'count'],
    index=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
).astype(float)


ROSTLAB_PROT_BERT = "Rostlab/prot_bert"


class DeconvolutionMethod(Enum):
    CONSTRAINED_EM = 'cem'
    EM = 'em'
    EMPIRICAL = 'empirical'
    LASSO = 'lasso'

    def __str__(self) -> str:
        return self.value


class DeconvolutionResult(Enum):
    CONFIDENT_HIT = 'confident_hit'
    CANDIDATE_HIT = 'candidate_hit'
    NOT_A_HIT = 'not_a_hit'

    def __str__(self) -> str:
        return self.value


class GenerateMode(Enum):
    GOLFY = 'golfy'
    CPSAT_SOLVER = 'cpsat_solver'

    def __str__(self) -> str:
        return self.value


class GolfyStrategy(Enum):
    GREEDY = 'greedy'
    RANDOM = 'random'
    REPEAT = 'repeat'

    def __str__(self) -> str:
        return self.value


class NumPlateWells(IntEnum):
    WELLS_24 = 24
    WELLS_48 = 48
    WELLS_96 = 96
    WELLS_384 = 384


class ReadoutFileType(Enum):
    POOL_IDS = 'pool_id'
    AID_PLATE_READER = 'aid_plate_reader'

    def __str__(self) -> str:
        return self.value


class SequenceSimilarityFunction(Enum):
    EUCLIDEAN = 'euclidean'
    COSINE = 'cosine'
    LEVENSHTEIN = 'levenshtein'

    def __str__(self) -> str:
        return self.value

