import solutions as sols
import numpy as np
import re
import itertools

from tqdm import tqdm

class Fasta_dna:
    def __init__(self, lines):
        assert lines[0][0] == '>'
        self.id = lines[0][1:].strip()
        self.dna_string = ""
        for line in lines[1:]:
            self.dna_string += line.strip()

    def __repr__(self):
        return self.id + " " + self.dna_string
    
    def count_acgt(self):
        acgt = [0, 0, 0, 0]
        for char in self.dna_string:
            if char in {"A", 'a'}:
                acgt[0] += 1
            elif char in {'C', 'c'}:
                acgt[1] += 1
            elif char in {'G', 'g'}:
                acgt[2] += 1
            elif char in {'T', 't'}:
                acgt[3] += 1
        return acgt[0], acgt[1], acgt[2], acgt[3]
    
    def to_rna(self):
        return self.dna_string.replace('T', 'U')
    
    def reverse_compliment(self):
        rev = self.dna_string[::-1]
        molecule_map = {"A":"T", "T":"A","C":"G","G":"C"}
        rev_comp = ""
        for char in rev:
            rev_comp += molecule_map[char]
        return rev_comp

    def gc_content(self):
        tot = len(self.dna_string)
        cgs = 0
        for char in self.dna_string:
            if char in ("C", 'G'):
                cgs += 1
        return (cgs/tot) * 100
    
    def motif_locations(self, substring):
        locations = []
        length = len(substring)
        for index in range(len(self.dna_string)):
            if self.dna_string[index:index+length] == substring:
                locations.append(index + 1)
        return locations
    
    def kmer_compostion(self, k):
        def string_from_list(lst):
            string = ""
            for char in lst:
                string += char
            return string
        alphabet = ["A", "C", "G", "T"]
        perms = list(itertools.product(alphabet, repeat=k))
        strs = [string_from_list(p) for p in perms]
        return [len(self.motif_locations(s)) for s in strs]
    
    def hamming_dist(self, comparison_string: str) -> int:
        if len(self.dna_string) != len(comparison_string):
            raise Exception("DNA lengths must be the same for hamming distance")
        else:
            diff = 0
            for index, char in enumerate(self.dna_string):
                if char != comparison_string[index]:
                    diff += 1
        return diff

    def translate_and_transcribe(self, introns):
        rna_string = self.to_rna()
        rna_introns = [intron.to_rna() for intron in introns]
        for rna_intron in rna_introns:
            rna_string = rna_string.replace(rna_intron, '', 1)
        return sols.codon_translation(rna_string)
    
    def reverse_palindromes(self):
        loc_lengths = []
        for i in range(len(self.dna_string)):
            for length in range(4, 13):
                substring = self.dna_string[i:i+length]
                if substring == sols.reverse_comp_string(substring) and i + length <= len(self.dna_string):
                    loc_lengths.append((i+1, length))
        return loc_lengths
    
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
                        self.dna_list.append(Fasta_dna(lines[start:end]))
                    start = index
            end = len(lines)
            self.dna_list.append(Fasta_dna(lines[start:end]))

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


