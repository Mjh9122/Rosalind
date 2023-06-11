import networkx as nx
import numpy as np
import math
import itertools
import Solutions as sols
import sys

from matplotlib import pyplot as plt
from tqdm import tqdm


class Fasta():
    def __init__(self, lines):
        """Initialize a Fasta object, consisting of a label, and a string (DNA, RNA, Amino Acids)

        Args:
            lines (List): First line is the label, starting with a >. The following lines are the string of the object
        """
        assert lines[0][0] == '>'
        self.id = lines[0][1:].strip()
        self.string = ""
        for line in lines[1:]:
            self.string += line.strip()
    
    def __repr__(self) -> str:
        """String representation, label and string

        Returns:
            str: see above
        """
        return self.id + " " + self.string
    
    def hamming_dist(self, comparison_string: str) -> int:
        """Calculate the hamming distance to a given string

        Args:
            comparison_string (str): The string to compare to

        Raises:
            Exception: Raised when the comparison string and self.string do not have the same length

        Returns:
            int: The number of characters that to not match when both strings are aligned by index
        """
        if len(self.string) != len(comparison_string):
            raise Exception("DNA lengths must be the same for hamming distance")
        else:
            diff = 0
            for index, char in enumerate(self.string):
                if char != comparison_string[index]:
                    diff += 1
        return diff
    
    def motif_locations(self, substring) -> list:
        """Finds the index of each appearance of the given substring

        Args:
            substring (str): substring to locate 

        Returns:
            list: list of indices where the substring appears
        """
        locations = []
        length = len(substring)
        for index in range(len(self.string)):
            if self.string[index:index+length] == substring:
                locations.append(index + 1)
        return locations
    
    def subsequence_indices(self, sequence_ob)->list:
        """Finds the indices of the subsequence in the Fasta object's string

        Args:
            sequence_ob (Fasta): The object whose string will be looked for as a subsequence

        Returns:
            list: list if indices that form the subsequence, if not all are found, returns list of as many as exist in the Fasta string
        """
        indices = []
        cur = 0
        for i, char in enumerate(self.string):
            try:
                if char == sequence_ob.string[cur]:
                    cur += 1
                    indices.append(i)
            except:
                pass
        return indices

    def common_subsequence(self, comparison_ob)->str:
        """Finds the longest common subsequence with another Fasta object

        Args:
            comparison_ob (Fasta): Fasta object whose string will be used to find the longest common subsequence

        Returns:
            str: the longest common subsequence (only one even if multiple are possible)
        """
        rows = self.string
        cols = comparison_ob.string
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

    def edit_distance(self, comparison_ob):
        """Finds the minimum edit distance between this and one other Fasta Object. 
        In this count, insertions, deletions, and swaps all count as 1 

        Args:
            comparison_ob (Fasta): Object whose string will be compared to

        Returns:
            int : number of operations needed to change Fasta string into comparison object's string
        """
        rows = self.string
        cols = comparison_ob.string
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
    
    def edit_dist_alignment(self, comparison_ob):
        """Returns the minimum edit distance and alignment of two strings
        Ex. 
        in: PRETTY PRTTEIN
        out: (4, PRETTY--, PR-TTEIN)

        Args:
            comparison_ob (Fasta): Object whose string will be compared to

        Returns:
            tuple: (edit distance, self.string with - included, comparison string with - included)
        """
        rows = self.string
        cols = comparison_ob.string
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
    
    def optimal_local_alignment(self, comparison_ob, mapper, gap_pen):
        """Finds the optimal local alignment of two strings, as defined by the mapping passed into the function

        Args:
            comparison_ob (Fasta): Object whose string will be compared to
            mapper (dict): dict whose keys should be the substitution chars, and values should be the score associated with them
            gap_pen (int): the penalty for insertions or deletions

        Returns:
            tuple: (max score, self.string substring associated with max score, comparison object substring associated with max score)
        """
        sys.setrecursionlimit(10**6)
        rows = self.string
        cols = comparison_ob.string
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
            """Recursively find the start of the alignment associated with the maximum alignment score

            Args:
                row (int): row
                col (int): column

            Returns:
                set: set of roots that could be the start of the max alignment
            """
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

        # Should never be reached (alignment not found that matches max score)
        return max_score, "failed"
    
    def fitting_alignment(self, motif_ob, mapper, gap_pen):
        """Finds the maximum alignment score of a motif in the Fasta string

        Args:
            motif_ob (Fasta): object whose string represents the motif
            mapper (dict): dict whose keys should be the substitution chars, and values should be the score associated with them
            gap_pen (int): the penalty for insertions or deletions

        Returns:
            int: maximum alinement score of the motif 
        """
        sys.setrecursionlimit(10**6)
        rows = motif_ob.string
        cols = self.string
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
        # modified edit dist alignment 
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