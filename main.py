import networkx as nx
import numpy as np
import math
import itertools

from matplotlib import pyplot as plt
from tqdm import tqdm


from Fasta_lib.Rna import Rna

with open('Inputs/rosalind_mmch(1).txt') as f:
    lines = f.readlines()
    rna_string = Rna(lines)
    ACGU = rna_string.count_acgu()
    print(ACGU)
    print(math.perm(max(ACGU[0], ACGU[3]), min(ACGU[0], ACGU[3])) * math.perm(max(ACGU[1], ACGU[2]), min(ACGU[1], ACGU[2])))
