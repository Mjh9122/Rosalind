import networkx as nx
import numpy as np
import math
import itertools
import Solutions as sols
import re

from Fasta import Fasta

from matplotlib import pyplot as plt
from tqdm import tqdm


class Dna(Fasta):
    """Class representing a DNA string in Fasta format

    Args:
        Fasta (_type_): parent class
    """
    def __init__(self, lines):
        """Init a Fasta DNA object

        Args:
            lines (List): First line is the label, starting with a >. The following lines are the string of the object
        """
        super().__init__(lines)

    def count_acgt(self):
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
            elif char in {'T', 't'}:
                acgt[3] += 1
        return acgt[0], acgt[1], acgt[2], acgt[3]
    
    def to_rna(self) -> str:
        """Turn the DNA into RNA using magic (turn T in U)

        Returns:
            str: RNA string
        """
        return self.string.replace('T', 'U')
    
    def reverse_compliment(self) -> str:
        """The reverse compliment of the DNA string
        reverse string and replace each letter with its base pair partner

        Returns:
            str: the reverse compliment string
        """
        rev = self.string[::-1]
        molecule_map = {"A":"T", "T":"A","C":"G","G":"C"}
        rev_comp = ""
        for char in rev:
            rev_comp += molecule_map[char]
        return rev_comp
    
    def gc_content(self) -> float:
        """Finds percent of DNA string that are either C or G

        Returns:
            float: Percent of DNA string that is C or G
        """
        tot = len(self.string)
        cgs = 0
        for char in self.string:
            if char in ("C", 'G'):
                cgs += 1
        return (cgs/tot) * 100
    
    def kmer_compostion(self, k) -> list:
        """Finds number of each k-mer in the DNA string

        Args:
            k (int): length of k-mer to find

        Returns:
            list: number of each k-mer (ordered lexicographically)
        """
        def string_from_list(lst):
            """Combines a list of strs into one str

            Args:
                lst (list): list of strings

            Returns:
                str: string of combined strings
            """
            string = ""
            for char in lst:
                string += char
            return string
        alphabet = ["A", "C", "G", "T"]
        perms = list(itertools.product(alphabet, repeat=k))
        strs = [string_from_list(p) for p in perms]
        return [len(self.motif_locations(s)) for s in strs]
    
    def translate_and_transcribe(self, introns):
        """Translate a DNA string into a protein string without the specific introns specified

        Args:
            introns (list): list of DNA substrings that should be ignored during the transcription process

        Returns:
            str: the protein str
        """
        rna_string = self.to_rna()
        rna_introns = [intron.to_rna() for intron in introns]
        for rna_intron in rna_introns:
            rna_string = rna_string.replace(rna_intron, '', 1)
        return sols.codon_translation(rna_string)
    
    def reverse_palindromes(self, min_length, max_length): 
        """Find reverse palindromes of specific lengths within the DNA string

        Args:
            min_length (int): minimum length
            max_length (int): maximum length

        Returns:
            list: list of tuples each tuple is the  length of the palindrome and the starting index of the palindrome
        """
        loc_lengths = []
        for i in range(len(self.string)):
            for length in range(min_length, max_length+1):
                substring = self.string[i:i+length]
                if substring == sols.reverse_comp_string(substring) and i + length <= len(self.string):
                    loc_lengths.append((i+1, length))
        return loc_lengths
    
    def transition_transversion_ratio(self, comparison_ob):
        """Find the ratio of transition vs translation mutations between another DNA string

        Args:
            comparison_ob (Fasta): Fasta object to compare to

        Returns:
            float: ratio of transition vs translation mutations
        """
        source = self.string
        target = comparison_ob.string
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
        """Finds all open reading frames in a DNA sequence and its reverse compliment by translating to RNA then a protein string, with each of the six possible start positions

        Returns:
            set: a set of valid proteins that are contained in the DNA string or its reverse compliment
        """
        reverse_comp = self.reverse_compliment()
        strs = [self.string, self.string[1:], self.string[2:], reverse_comp, reverse_comp[1:], reverse_comp[2:]]
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
    