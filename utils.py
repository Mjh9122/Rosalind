import networkx as nx
import numpy as np
import math
import itertools

from Dna import Dna
from Rna import Rna
from Protein import Protein
from Fasta import Fasta
import Solutions as sols

from matplotlib import pyplot as plt
from tqdm import tqdm


def blossom_map():
    with open('Solutions/utils/blossum62.txt') as f:
        lines = f.readlines()
        matrix = [line.strip().split() for line in lines]
        top_labels = matrix[0]
        matrix = matrix[1:]
        side_labels = [matrix[i][0] for i in range(len(matrix))]        
        matrix = [matrix[i][1:] for i in range(len(matrix))]
        blossom_map = dict()
        for i, top in enumerate(top_labels):
            for j, side in enumerate(side_labels):
                blossom_map[top+side] = int(matrix[j][i])

        return blossom_map
    
def pam_map():
    with open('Solutions/utils/pam50.txt') as f:
        lines = f.readlines()
        matrix = [line.strip().split() for line in lines]
        top_labels = matrix[0]
        matrix = matrix[1:]
        side_labels = [matrix[i][0] for i in range(len(matrix))]        
        matrix = [matrix[i][1:] for i in range(len(matrix))]
        blossom_map = dict()
        for i, top in enumerate(top_labels):
            for j, side in enumerate(side_labels):
                blossom_map[top+side] = int(matrix[j][i])

        return blossom_map
    
def mismatch_map(match_score, missmatch_score):
    with open('Solutions/utils/blossum62.txt') as f:
        lines = f.readlines()
        matrix = [line.strip().split() for line in lines]
        top_labels = matrix[0]
        side_labels = [matrix[i][0] for i in range(len(matrix))]        
        mismatch_map = dict()
        for i, top in enumerate(top_labels):
            for j, side in enumerate(side_labels):
                if top == side:
                    mismatch_map[top+side] = 1
                else:
                    mismatch_map[top+side] = -1
        return mismatch_map
    
int_to_nucleotide = {0:'A', 1:"C", 2:"G", 3:'T'}
nucleotide_to_int = {'A':0, "C":1, "G":2, 'T':3}