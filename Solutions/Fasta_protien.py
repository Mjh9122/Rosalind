import networkx as nx
import numpy as np
import math
import itertools
import Fasta
import Solutions as sols

from matplotlib import pyplot as plt
from tqdm import tqdm


class Fasta_protein(Fasta):
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
