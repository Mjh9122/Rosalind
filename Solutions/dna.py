import solutions as sols
import numpy as np
import re
import itertools
import sys

from tqdm import tqdm

class Fasta:
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
    
    def subsequence_indices(self, sequence_ob):
        indices = []
        cur = 0
        for i, char in enumerate(self.dna_string):
            try:
                if char == sequence_ob.dna_string[cur]:
                    cur += 1
                    indices.append(i)
            except:
                pass
        return indices
    
    def common_subsequence(self, comparison_ob):
        rows = self.dna_string
        cols = comparison_ob.dna_string
        longest_substr = [["" for _ in range(len(cols)+1)] for _ in range(len(rows)+1)]
        for r in range(1, len(rows)+1):
            for c in range(1, len(cols)+1):
                if rows[r-1] == cols[c-1]:
                    longest_substr[r][c] = longest_substr[r-1][c-1] + rows[r-1]
                elif len(longest_substr[r][c-1]) > len(longest_substr[r-1][c]):
                    longest_substr[r][c] = longest_substr[r][c-1]
                else: 
                    longest_substr[r][c] = longest_substr[r-1][c]
        return longest_substr[-1][-1]

    def transition_transversion_ratio(self, comparison):
        source = self.dna_string
        target = comparison.dna_string
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

    def edit_distance(self, comparison_ob):
        rows = self.dna_string
        cols = comparison_ob.dna_string
        dist = [[0 for _ in range(len(cols)+1)] for _ in range(len(rows)+1)]
        for i in range(len(rows) + 1):
            dist[i][0] = i
        for j in range(len(cols)+1):
            dist[0][j] = j
        for r in range(1, len(rows)+1):
            for c in range(1, len(cols)+1):
                if rows[r-1] == cols[c-1]:
                    dist[r][c] = dist[r-1][c-1]
                else:
                    dist[r][c] = min(dist[r-1][c], dist[r-1][c-1], dist[r][c-1]) + 1
        return dist[-1][-1]
    
    def edit_dist_allignment(self, comparison_ob):
        rows = self.dna_string
        cols = comparison_ob.dna_string
        # (ROW, COL, NUM)
        dist = [[("", "", 0) for _ in range(len(cols)+1)] for _ in range(len(rows)+1)]
        for i in range(1, len(rows) + 1):
            dist[i][0] = (dist[i-1][0][0]+rows[i-1], dist[i-1][0][1]+"-", i)
        for j in range(1, len(cols)+1):
            dist[0][j] = (dist[0][j-1][0]+"-", dist[0][j-1][1]+cols[j-1], j)
        for r in range(1, len(rows)+1):
            for c in range(1, len(cols)+1):
                if rows[r-1] == cols[c-1]:
                    dist[r][c] = (dist[r-1][c-1][0] + rows[r-1], dist[r-1][c-1][1] + rows[r-1], dist[r-1][c-1][2])
                elif dist[r-1][c-1][2] <= dist[r-1][c][2] and dist[r-1][c-1][2] <= dist[r][c-1][2]:
                    dist[r][c] = (dist[r-1][c-1][0] + rows[r-1], dist[r-1][c-1][1] + cols[c-1], dist[r-1][c-1][2]+1)
                elif dist[r-1][c][2] <= dist[r][c-1][2]:
                    if '-' in dist[r-1][c][0] and dist[r-1][c][0][-1] == '-':
                        dist[r][c] = (dist[r-1][c][0][:-1] + rows[r-1], dist[r-1][c][1], dist[r-1][c][2])
                    else:
                        dist[r][c] = (dist[r-1][c][0] + rows[r-1], dist[r-1][c][1] + "-", dist[r-1][c][2]+1)
                else:
                    if '-' in dist[r][c-1][0] and dist[r][c-1][1][-1] == '-':
                        dist[r][c] = (dist[r][c-1][0], dist[r][c-1][1][:-1]+cols[c-1], dist[r][c-1][2])
                    else:
                        dist[r][c] = (dist[r][c-1][0] + "-", dist[r][c-1][1] + cols[c-1], dist[r][c-1][2]+1)
        return dist[-1][-1]
    
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
    
    def optimal_local_alignment(self, comparison_ob, mapper, gap_pen):
        sys.setrecursionlimit(10**6)
        rows = self.dna_string
        cols = comparison_ob.dna_string
        dist = [[0 for _ in range(len(cols)+1)] for _ in range(len(rows)+1)]
        traceback = [[[(0, 0)] for _ in range(len(cols) + 1)] for _ in range(len(rows)+1)]
        for r in range(1, len(rows)+1):
            for c in range(1, len(cols)+1):
                next_cell_val = max(dist[r-1][c-1] + mapper[rows[r-1]+cols[c-1]], dist[r-1][c] - gap_pen, dist[r][c-1] - gap_pen, 0)
                dist[r][c] = next_cell_val
                traceback[r][c] = []
                if dist[r-1][c-1] + mapper[rows[r-1]+cols[c-1]] == next_cell_val:
                    traceback[r][c].append((-1, -1))
                if dist[r-1][c] - gap_pen == next_cell_val:
                    traceback[r][c].append((-1, 0))
                if dist[r][c-1] - gap_pen == next_cell_val:
                    traceback[r][c].append((0, -1))
                if 0 == next_cell_val:
                    traceback[r][c].append((0, 0))

        max_score = 0
        max_loc = (0, 0)
        for row_index, row in enumerate(dist):
            if max(row) > max_score:
                max_row = row_index
                max_col = row.index(max(row))
                max_loc = max_row, max_col
                max_score = max(row)

        end_row, end_col = max_loc

        def find_roots(row, col):
            roots = set()
            if dist[row][col] == 0:
                roots.add((row, col))
                return roots
            else:
                for direction in traceback[row][col]:
                    return roots.union(find_roots(row + direction[0], col + direction[1]))

        starts = find_roots(end_row, end_col)
        for (start_row, start_col) in list(starts):
            return max_score, rows[start_row:end_row], cols[start_col:end_col]

        return max_score, "failed"
    
    def fitting_alignment(self, motif_ob, mapper, gap_pen):
        sys.setrecursionlimit(10**6)
        rows = motif_ob.dna_string
        cols = self.dna_string
        dist = [[0 for _ in range(len(cols)+1)] for _ in range(len(rows)+1)]
        for i in range(len(rows)+1):
            dist[i][0] = i * -gap_pen
        traceback = [[[(0, 0)] for _ in range(len(cols) + 1)] for _ in range(len(rows)+1)]
        for r in range(1, len(rows)+1):
            for c in range(1, len(cols)+1):
                next_cell_val = max(dist[r-1][c-1] + mapper[rows[r-1]+cols[c-1]], dist[r-1][c] - gap_pen, dist[r][c-1] - gap_pen)
                dist[r][c] = next_cell_val
                traceback[r][c] = []
                if dist[r-1][c-1] + mapper[rows[r-1]+cols[c-1]] == next_cell_val:
                    traceback[r][c].append((-1, -1))
                if dist[r-1][c] - gap_pen == next_cell_val:
                    traceback[r][c].append((-1, 0))
                if dist[r][c-1] - gap_pen == next_cell_val:
                    traceback[r][c].append((0, -1))

        max_score = max(dist[-1])
        max_loc = len(rows), dist[-1].index(max(dist[-1]))

        end_row, end_col = max_loc

        def find_roots(row, col):
            roots = set()
            if row == 0:
                roots.add((row, col))
                return roots
            else:
                for direction in traceback[row][col]:
                    return roots.union(find_roots(row + direction[0], col + direction[1]))

        starts = find_roots(end_row, end_col)
        for (start_row, start_col) in list(starts):
            col_string = cols[start_col:end_col]

        ###########################################
        cols = col_string
        # (ROW, COL, NUM)
        dist = [[("", "", 0) for _ in range(len(cols)+1)] for _ in range(len(rows)+1)]
        for i in range(1, len(rows) + 1):
            dist[i][0] = (dist[i-1][0][0]+rows[i-1], dist[i-1][0][1]+"-", -i)
        for j in range(1, len(cols)+1):
            dist[0][j] = (dist[0][j-1][0]+"-", dist[0][j-1][1]+cols[j-1], -j)
        for r in range(1, len(rows)+1):
            for c in range(1, len(cols)+1):
                if rows[r-1] == cols[c-1]:
                    dist[r][c] = (dist[r-1][c-1][0] + rows[r-1], dist[r-1][c-1][1] + rows[r-1], dist[r-1][c-1][2]+1)
                elif dist[r-1][c-1][2] >= dist[r-1][c][2] and dist[r-1][c-1][2] >= dist[r][c-1][2]:
                    dist[r][c] = (dist[r-1][c-1][0] + rows[r-1], dist[r-1][c-1][1] + cols[c-1], dist[r-1][c-1][2]-1)
                elif dist[r-1][c][2] >= dist[r][c-1][2]:
                    if '-' in dist[r-1][c][0] and dist[r-1][c][0][-1] == '-':
                        dist[r][c] = (dist[r-1][c][0][:-1] + rows[r-1], dist[r-1][c][1], dist[r-1][c][2])
                    else:
                        dist[r][c] = (dist[r-1][c][0] + rows[r-1], dist[r-1][c][1] + "-", dist[r-1][c][2]-1)
                else:
                    if '-' in dist[r][c-1][0] and dist[r][c-1][1][-1] == '-':
                        dist[r][c] = (dist[r][c-1][0], dist[r][c-1][1][:-1]+cols[c-1], dist[r][c-1][2])
                    else:
                        dist[r][c] = (dist[r][c-1][0] + "-", dist[r][c-1][1] + cols[c-1], dist[r][c-1][2]-1)
        ##########################################

        return dist[-1][-1]
        
    


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


