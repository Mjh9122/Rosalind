def dna_count_acgt(dna_string):
    acgt = [0, 0, 0, 0]
    for char in dna_string:
        if char in {"A", 'a'}:
            acgt[0] += 1
        elif char in {'C', 'c'}:
            acgt[1] += 1
        elif char in {'G', 'g'}:
            acgt[2] += 1
        elif char in {'T', 't'}:
            acgt[3] += 1
    return acgt[0], acgt[1], acgt[2], acgt[3]
        
def dna_to_rna(dna_string):
    return dna_string.replace('T', 'U')

def dna_reverse_compliment(dna_string):
    rev = dna_string[::-1]
    molecule_map = {"A":"T", "T":"A","C":"G","G":"C"}
    rev_comp = ""
    for char in rev:
        rev_comp += molecule_map[char]
    return rev_comp
