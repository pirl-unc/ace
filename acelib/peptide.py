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
The purpose of this python3 script is to set up a Peptide class in ACE.
This will assist in the calculation of sequence level features to be included in CP_SAT.
"""


import re
from sklearn.metrics.pairwise import cosine_similarity


class Peptide:

    def __init__(self, id, sequence, tokenizer, model):
        self.id = id
        self.sequence = re.sub(r'[UZOB]', 'X', sequence)

    def __str__(self):
        return self.sequence

    def embed(self, sequence, tokenizer, encoder):
        sequence_w_spaces = str.replace(sequence,  "", " ")[1:-1]
        encoded_input = tokenizer(sequence_w_spaces, return_tensors='pt')
        output = encoder(**encoded_input)
        return output['pooler_output'].detach().numpy()

    def distance2(self, peptide2):
        sim = cosine_similarity(self.embedding, peptide2.embedding).item(0)
        return sim
