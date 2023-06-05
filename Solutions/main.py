import solutions as sols
import dna
import numpy as np
import itertools

from tqdm import tqdm


lst = dna.Fasta_List_Ops('Inputs/rosalind_tran.txt').dna_list
print(lst[0].transition_transversion_ratio(lst[1]))
    