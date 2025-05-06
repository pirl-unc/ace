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
The purpose of this python3 script is to implement the HitPeptide class.
"""


from dataclasses import dataclass, field
from typing import List, Tuple
from .constants import DeconvolutionResult
from .logger import get_logger
from .peptide import Peptide


logger = get_logger(__name__)


@dataclass
class DeconvolvedPeptide:
    peptide: Peptide
    estimated_spot_count: float
    result: DeconvolutionResult
    hit_pool_ids: List[int] = field(default_factory=list)

    @property
    def id(self) -> str:
        return self.peptide.id

    @property
    def sequence(self) -> str:
        return self.peptide.sequence
