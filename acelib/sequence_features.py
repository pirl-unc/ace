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
Purpose of this Code is to Handle the Sequence Similarity Elements of the ACE CP-SAT
"""


import torch
import numpy as np
import difflib


def bert_embed(sequence, tokenizer, encoder):
    """

    Parameters
    ----------
    sequence    :
    tokenizer   :
    encoder     :

    Returns
    -------
    """
    encoded_input = tokenizer(" ".join(sequence), return_tensors='pt')
    output = encoder(**encoded_input)
    return torch.mean(output.last_hidden_state, 1).detach().numpy()[0]


def cosine_similarity(emb1, emb2):
    return np.dot(emb1, emb2) / (np.linalg.norm(emb1)*np.linalg.norm(emb2))


def dummy_embed(sequence):
    return sequence


def kmer_similarity(seq1, seq2):
    return difflib.SequenceMatcher(None, seq1, seq2).find_longest_match().size


def ratio_similarity(seq1, seq2):
    return difflib.SequenceMatcher(None, seq1, seq2).ratio()


def is_disallowed(seq1, seq2, emb_fxn, sim_fxn, threshold):
    if sim_fxn(emb_fxn(seq1), emb_fxn(seq2)) >= threshold:
        return True
    else:
        return False


def get_disallowed_peptides(seqs, emb_fxn, sim_fxn, threshold):
    return [(idx, idx2) for idx, seq1 in enumerate(seqs) for idx2, seq2 in enumerate(seqs[idx + 1:]) if is_disallowed(seq1, seq2, emb_fxn, sim_fxn, threshold)]
