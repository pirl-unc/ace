"""
The purpose of this python3 script is to set up a Peptide class in ACE.
This will assist in the calculation of sequence level features to be included in CP_SAT.

Last updated date: October 12, 2022

Author: Jin Seok (Andy) Lee, Dhuvarakesh Karthikeyan
"""
import re
from transformers import BertModel, BertTokenizer


class Peptide:

    def __init__(self, id, sequence, tokenizer, model):
        self.id = id
        self.amino_acids = list('GAPVLIMYFWSTCNQKRHDEX')
        self.sequence = re.sub(r'[UZOB]', 'X', sequence)
        self.aa2index = {}
        self.aa2count = dict(zip(self.amino_acids, [0]*21))
        self.index2aa = {0:'*', 1:'#'} #Start Index and Stop Index
        self.n_residues = 2
        self.__load_sequence__(self.sequence)
        self.embedding = self.embed(self.sequence, tokenizer, model)

    def __load_sequence__(self, processed_sequence):
        for aa in list(processed_sequence):
            self.__add_residue__(aa)

    def __add_residue__(self, aa):
        if aa not in self.aa2index:
            self.aa2index[aa] = self.n_residues
            self.index2aa[self.n_residues] = aa
            self.n_residues += 1

        self.aa2count[aa] += 1

    def __str__(self):
        return self.sequence

    def embed(self, sequence, tokenizer, encoder):
        encoded_input = tokenizer(sequence, return_tensors='pt')
        output = encoder(**encoded_input)
        return output
