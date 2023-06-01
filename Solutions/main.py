import solutions as sols
import dna
import numpy as np
import itertools

lst = dna.Fasta_List_Ops('Inputs/tests.txt').dna_list
for item in lst[0].optimal_local_alignment(lst[1], sols.pam_map, 5):
    print(item)