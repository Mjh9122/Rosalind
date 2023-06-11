import networkx as nx
import numpy as np
import math
import itertools
from Fasta import Fasta
import Solutions as sols

from matplotlib import pyplot as plt
from tqdm import tqdm


amino_acid_weight_map = {
    "A":71.03711,
    "C":103.00919,
    "D":115.02694,
    "E":129.04259,
    "F":147.06841,
    "G":57.02146,
    "H":137.05891,
    "I":113.08406,
    "K":128.09496,
    "L":113.08406,
    "M":131.04049,
    "N":114.04293,
    "P":97.05276,
    "Q":128.05858,
    "R":156.10111,
    "S":87.03203,
    "T":101.04768,
    "V":99.06841,
    "W":186.07931,
    "Y":163.06333,
}

class Protein(Fasta):
    def __init__(self, lines):
        super().__init__(lines)

    # Returns the minimum weighted edit distance, where the weights for each swap is found in mapper
    def edit_distance_weighted(self, comparison_ob, mapper):
        rows = self.dna_string
        cols = comparison_ob.dna_string
        dist = [[0 for _ in range(len(cols)+1)] for _ in range(len(rows)+1)]
        for i in range(len(rows) + 1):
            dist[i][0] = i * -5
        for j in range(len(cols)+1):
            dist[0][j] = j * -5
        for r in range(1, len(rows)+1):
            for c in range(1, len(cols)+1):
                if rows[r-1] == cols[c-1]:
                    dist[r][c] = dist[r-1][c-1] + mapper[rows[r-1]+cols[c-1]]
                else:
                    dist[r][c] = max(dist[r-1][c] - 5, dist[r-1][c-1] + mapper[rows[r-1]+cols[c-1]], dist[r][c-1] - 5)
        return dist[-1][-1]
    
    def protein_mass(self):
        mass = 0
        for char in self.string:
            mass += amino_acid_weight_map[char]
        return mass
