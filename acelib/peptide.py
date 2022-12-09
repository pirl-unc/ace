"""
The purpose of this python3 script is to set up a Peptide class in ACE.
This will assist in the calculation of sequence level features to be included in CP_SAT.

Last updated date: October 12, 2022

Author: Jin Seok (Andy) Lee, Dhuvarakesh Karthikeyan
"""
import re
import numpy as np
import difflib

def bert_embed(self, sequence, tokenizer, encoder):
    encoded_input = tokenizer(sequence, return_tensors='pt')
    output = encoder(**encoded_input)
    return output


def cosine_similarity(emb1, emb2):
    return np.dot(emb1, emb2) / (np.linalg.norm(emb1)*np.linalg.norm(emb2))


def dummy_embed(sequence):
    return sequence


def kmer_similarity(seq1, seq2):
    return difflib.SequenceMatcher(None, seq1, seq2).find_longest_match().size


def ratio_similarity(seq1, seq2):
    return difflib.SequenceMatcher(None, seq1, seq2).ratio()


class Peptide:

    def __init__(self, id, sequence, tokenizer, model):
        self.id = id
        self.sequence = re.sub(r'[UZOB]', 'X', sequence)

    def __str__(self):
        return self.sequence
