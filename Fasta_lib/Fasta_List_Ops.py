import Fasta
import networkx as nx
import numpy as np
import math
import itertools
import Solutions as sols

from matplotlib import pyplot as plt
from tqdm import tqdm


class Fasta_List_Ops:
    def __init__(self, filepath):
        self.dna_list = []
        with open(filepath) as f:
            lines = f.readlines()
            start, end = 0, 0
            for index, line in enumerate(lines):
                if line[0] == '>':
                    end = index
                    if start != end:
                        self.dna_list.append(Fasta(lines[start:end]))
                    start = index
            end = len(lines)
            self.dna_list.append(Fasta(lines[start:end]))

    def largest_substring(self):
        srted = sorted(self.dna_list, key=lambda x: len(x.dna_string))
        substrs = list()
        for i in range(len(srted[0].dna_string)):
            for j in range(i, len(srted[0].dna_string)):
                substring = srted[0].dna_string[i:j]
                in_all = True
                for dna_ob in self.dna_list:
                    if substring not in dna_ob.dna_string:
                        in_all = False
                        break
                if in_all:
                    substrs.append(substring)

        print(max(substrs, key=lambda x: len(x)))

    def profile_matrix(self):
        length = len(max(self.dna_list, key=lambda x: len(x.dna_string)).dna_string)
        profile_mat = np.zeros((4, length), dtype=np.int64)
        for dna_ob in self.dna_list:
            for index, char in enumerate(dna_ob.dna_string):
                if char == 'A':
                    profile_mat[0][index] += 1
                elif char == 'C':
                    profile_mat[1][index] += 1
                elif char == 'G':
                    profile_mat[2][index] += 1
                elif char == 'T':
                    profile_mat[3][index] += 1
        return profile_mat

    def consensus_string(self):
        profile_matrix = self.profile_matrix()
        max_indices = [(np.argmax(profile_matrix[:,j])) for j in range(len(profile_matrix[0]))]
        consensus_string = ""
        for i in max_indices:
            consensus_string += sols.int_to_nucleotide[i]
        return consensus_string
        

    def overlap_graph(self, overlap_len):
        edges = []
        for source in self.dna_list:
            for target in self.dna_list:
                if source != target:
                    if source.dna_string[-overlap_len:] == target.dna_string[:overlap_len]:
                        edges.append((source.id, target.id))
        return edges
    
    def __repr__(self) -> str:
        return str([dna.id for dna in self.dna_list])


