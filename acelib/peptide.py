"""
The purpose of this python3 script is to set up a Peptide class in ACE.
This will assist in the calculation of sequence level features to be included in CP_SAT.

Last updated date: October 12, 2022

Author: Jin Seok (Andy) Lee, Dhuvarakesh Karthikeyan
"""
import re

class Peptide:

    def __init__(self, id, sequence, embedding):
        self.id = id
        self.amino_acids = split(['GAPVLIMYFWSTCNQKRHDEX'])
        self.aa2index = {}
        self.aa2count = dict(zip(self.amino_acids, [0]*21))
        self.index2aa = {0:'*', 1:'#'} #Start Index and Stop Index
        self.n_residues = 2

    def load_sequence(self, sequence):
        processed_sequence = re.sub(r'[UZOB]', 'X', sequence)
        for aa in split(processed_sequence):
            self.add_residue(aa)

    def add_residue(self, aa):
        if aa not in self.aa2index:
            self.aa2index[aa] = self.n_residues
            self.index2aa[self.n_residues] = aa
            self.n_residues += 1
        else:
            continue
        self.aa2count[aa] += 1
