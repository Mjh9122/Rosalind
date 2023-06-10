import Fasta
import networkx as nx
import numpy as np
import math
import itertools
import solutions as sols
import re

from matplotlib import pyplot as plt
from tqdm import tqdm


class Fasta_dna(Fasta):
    def __init__(self, lines):
        super().__init__(lines)

    def count_acgt(self):
        acgt = [0, 0, 0, 0]
        for char in self.string:
            if char in {"A", 'a'}:
                acgt[0] += 1
            elif char in {'C', 'c'}:
                acgt[1] += 1
            elif char in {'G', 'g'}:
                acgt[2] += 1
            elif char in {'T', 't'}:
                acgt[3] += 1
        return acgt[0], acgt[1], acgt[2], acgt[3]
    
    def to_rna(self) -> str:
        return self.string.replace('T', 'U')
    
    def reverse_compliment(self) -> str:
        rev = self.string[::-1]
        molecule_map = {"A":"T", "T":"A","C":"G","G":"C"}
        rev_comp = ""
        for char in rev:
            rev_comp += molecule_map[char]
        return rev_comp
    
    def gc_content(self) -> float:
        tot = len(self.string)
        cgs = 0
        for char in self.string:
            if char in ("C", 'G'):
                cgs += 1
        return (cgs/tot) * 100
    
    def kmer_compostion(self, k) -> list:
        def string_from_list(lst):
            string = ""
            for char in lst:
                string += char
            return string
        alphabet = ["A", "C", "G", "T"]
        perms = list(itertools.product(alphabet, repeat=k))
        strs = [string_from_list(p) for p in perms]
        return [len(self.motif_locations(s)) for s in strs]
    
    def translate_and_transcribe(self, introns):
        rna_string = self.to_rna()
        rna_introns = [intron.to_rna() for intron in introns]
        for rna_intron in rna_introns:
            rna_string = rna_string.replace(rna_intron, '', 1)
        return sols.codon_translation(rna_string)
    
    # List containing the location and length of every reverse palindrome between a min and max length 
    def reverse_palindromes(self, min_length, max_length):
        loc_lengths = []
        for i in range(len(self.string)):
            for length in range(min_length, max_length+1):
                substring = self.string[i:i+length]
                if substring == sols.reverse_comp_string(substring) and i + length <= len(self.string):
                    loc_lengths.append((i+1, length))
        return loc_lengths
    
    def transition_transversion_ratio(self, comparison):
        source = self.string
        target = comparison.string
        assert len(source) == len(target)
        transitions = 0
        transversions = 0
        for i, char1 in enumerate(source):
            if char1 != (char2 := target[i]):
                if char1 in ('A', "G") and char2 in ('A', 'G') or \
                char1 in ('C', 'T') and char2 in ('C', 'T'):
                    transitions += 1
                else:
                    transversions += 1
        return transitions/transversions
    
    def open_reading_frames(self):
        reverse_comp = self.reverse_compliment()
        strs = [self.dna_string, self.dna_string[1:], self.dna_string[2:], reverse_comp, reverse_comp[1:], reverse_comp[2:]]
        rna_strs = [sols.to_rna(s) for s in strs]
        codons_strs = [sols.codon_translation(rna) for rna in rna_strs]
        valid_start_strings = []
        for codon_str in codons_strs:
            try:
                start_indeces = [m.start() for m in re.finditer("M", codon_str)]
                for i in start_indeces:
                    valid_start_strings.append(codon_str[i:])
            except:
                continue
        valid_proteins = [cod[:cod.index('Stop')] for cod in valid_start_strings if 'Stop' in cod]
        return set(valid_proteins)
    