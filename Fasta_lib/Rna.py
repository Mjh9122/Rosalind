from Fasta_lib.Fasta import Fasta
import networkx as nx
import numpy as np
import math
import itertools
import Solutions as sols


from matplotlib import pyplot as plt
from tqdm import tqdm

codon_map = {
        "UUU": "F",
        "CUU": "L",
        "AUU": "I",
        "GUU": "V",
        "UUC": "F",
        "CUC": "L",
        "AUC": "I",
        "GUC": "V",
        "UUA": "L",
        "CUA": "L",
        "AUA": "I",
        "GUA": "V",
        "UUG": "L",
        "CUG": "L",
        "AUG": "M",
        "GUG": "V",
        "UCU": "S",
        "CCU": "P",
        "ACU": "T",
        "GCU": "A",
        "UCC": "S",
        "CCC": "P",
        "ACC": "T",
        "GCC": "A",
        "UCA": "S",
        "CCA": "P",
        "ACA": "T",
        "GCA": "A",
        "UCG": "S",
        "CCG": "P",
        "ACG": "T",
        "GCG": "A",
        "UAU": "Y",
        "CAU": "H",
        "AAU": "N",
        "GAU": "D",
        "UAC": "Y",
        "CAC": "H",
        "AAC": "N",
        "GAC": "D",
        "UAA": "Stop",
        "CAA": "Q",
        "AAA": "K",
        "GAA": "E",
        "UAG": "Stop",
        "CAG": "Q",
        "AAG": "K",
        "GAG": "E",
        "UGU": "C",
        "CGU": "R",
        "AGU": "S",
        "GGU": "G",
        "UGC": "C",
        "CGC": "R",
        "AGC": "S",
        "GGC": "G",
        "UGA": "Stop",
        "CGA": "R",
        "AGA": "R",
        "GGA": "G",
        "UGG": "W",
        "CGG": "R",
        "AGG": "R",
        "GGG": "G",
    }

class Rna(Fasta):
    def __init__(self, lines):
        super().__init__(lines)

    def count_acgu(self):
        """Counts the appearance of each letter in the DNA string

        Returns:
            tuple: (num of As, num of Cs, num of Gs, num of Ts)
        """
        acgt = [0, 0, 0, 0]
        for char in self.string:
            if char in {"A", 'a'}:
                acgt[0] += 1
            elif char in {'C', 'c'}:
                acgt[1] += 1
            elif char in {'G', 'g'}:
                acgt[2] += 1
            elif char in {'U', 'u'}:
                acgt[3] += 1
        return acgt[0], acgt[1], acgt[2], acgt[3]
    
    def to_protein_string(self):
        amino_acid_string = ""
        codons = [self.string[i:i+3] for i in range(0,len(self.string),3)]
        for codon in codons[:-1]:
            try:
                amino_acid_string += codon_map[codon]
            except:
                return ""
        return amino_acid_string
    
    def codon_combinations(self, modulus):
        reverse_codon_count = dict()
        for acid in codon_map.values():
            if acid in reverse_codon_count:
                reverse_codon_count[acid] += 1
            else:
                reverse_codon_count[acid] = 1
        total = 1
        for char in self.string:
            total *= reverse_codon_count[char]
        total *= reverse_codon_count['Stop']
        return total % 1_000_000