import solutions as sols
import dna

with open('Inputs/rosalind_mrna.txt') as f:
    string = f.readline().strip()
    print(sols.codon_combinations(string))