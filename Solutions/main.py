import solutions as sols
import dna
import numpy as np



lst = dna.Fasta_List_Ops('Inputs/rosalind_cons.txt')
print(lst.consensus_string())
matrix = lst.profile_matrix()

for i, row in enumerate(matrix):
    row_str = sols.int_to_nucleotide[i]+":"
    for num in row: 
        row_str += " " + str(num)
    print(row_str)

