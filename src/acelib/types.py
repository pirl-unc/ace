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
The purpose of this python3 script is to implement data types used in ACE.
"""


import pandas as pd
from typing import Iterable, Mapping, Tuple


CoverageId = int
DeconvolutionLabel = str
PeptideId = str
PeptideSequence = str
MaxPeptidesPerBlock = int
NumPeptides = int
NumPeptidesPerPool = int
NumCoverage = int
PeptideSequenceSimilarityScore = float
PlateId = int
PoolId = int
SpotCount = int
WellId = str

PeptidePairs = Iterable[Tuple[PeptideId, PeptideId]]
Peptides = Iterable[Tuple[PeptideId, PeptideSequence]]
Assignments = Mapping[CoverageId, Mapping[PoolId, Iterable[Tuple[PeptideId, PeptideSequence]]]]
HitPeptides = Iterable[Tuple[PeptideId, PeptideSequence, DeconvolutionLabel, Iterable[PoolId]]]
PlateReadoutObservations = Iterable[Tuple[PlateId, WellId, SpotCount]]

