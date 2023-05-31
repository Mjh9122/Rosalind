import numpy as np

int_to_nucleotide = {0:'A', 1:"C", 2:"G", 3:'T'}
nucleotide_to_int = {'A':0, "C":1, "G":2, 'T':3}


codon_map = {
        "UUU": "F",
        "CUU": "L",
        "AUU": "I",
        "GUU": "V",
        "UUC": "F",
        "CUC": "L",
        "AUC": "I",
        "GUC": "V",
        "UUA": "L",
        "CUA": "L",
        "AUA": "I",
        "GUA": "V",
        "UUG": "L",
        "CUG": "L",
        "AUG": "M",
        "GUG": "V",
        "UCU": "S",
        "CCU": "P",
        "ACU": "T",
        "GCU": "A",
        "UCC": "S",
        "CCC": "P",
        "ACC": "T",
        "GCC": "A",
        "UCA": "S",
        "CCA": "P",
        "ACA": "T",
        "GCA": "A",
        "UCG": "S",
        "CCG": "P",
        "ACG": "T",
        "GCG": "A",
        "UAU": "Y",
        "CAU": "H",
        "AAU": "N",
        "GAU": "D",
        "UAC": "Y",
        "CAC": "H",
        "AAC": "N",
        "GAC": "D",
        "UAA": "Stop",
        "CAA": "Q",
        "AAA": "K",
        "GAA": "E",
        "UAG": "Stop",
        "CAG": "Q",
        "AAG": "K",
        "GAG": "E",
        "UGU": "C",
        "CGU": "R",
        "AGU": "S",
        "GGU": "G",
        "UGC": "C",
        "CGC": "R",
        "AGC": "S",
        "GGC": "G",
        "UGA": "Stop",
        "CGA": "R",
        "AGA": "R",
        "GGA": "G",
        "UGG": "W",
        "CGG": "R",
        "AGG": "R",
        "GGG": "G",
    }

amino_acid_weight_map = {
    "A":71.03711,
    "C":103.00919,
    "D":115.02694,
    "E":129.04259,
    "F":147.06841,
    "G":57.02146,
    "H":137.05891,
    "I":113.08406,
    "K":128.09496,
    "L":113.08406,
    "M":131.04049,
    "N":114.04293,
    "P":97.05276,
    "Q":128.05858,
    "R":156.10111,
    "S":87.03203,
    "T":101.04768,
    "V":99.06841,
    "W":186.07931,
    "Y":163.06333,
}

def immortal_rabbit_pop(months, rabbits_per_litter):
    pop, pop_past = 1, 1
    for _ in np.arange(months - 2):
        pop, pop_past = pop_past * rabbits_per_litter + pop, pop
    return pop

def mortal_rabbit_pop(months, life_span):
    pop = [0 for _ in range(life_span)]
    pop[-1] = 1
    for month in range(2, months+1):
        pop.append(sum(pop[:-1]))
        pop = pop[1:]
    return sum(pop)

def dominant_allele_by_population(homo_pos, het, homo_neg):
    tot = homo_neg + homo_pos + het
    pos_first_parent = homo_pos / tot
    het_first_parent = het / tot
    neg_first_parent = homo_neg / tot

    # Pos + X
    pos_chance = pos_first_parent

    # Het + X
    pos_chance += het_first_parent * (homo_pos / (tot - 1))
    pos_chance += het_first_parent * ((het - 1) / (tot - 1)) * 0.75
    pos_chance += het_first_parent * (homo_neg / (tot - 1)) * 0.5

    # Neg + X
    pos_chance += neg_first_parent * (homo_pos / (tot - 1))
    pos_chance += neg_first_parent * ((het) / (tot - 1)) * 0.5

    return pos_chance

def dominant_allele_by_parrent_matching(pos_pos, pos_het, pos_neg, het_het, het_neg, neg_neg):
    return pos_pos * 2 + pos_het * 2 + pos_neg * 2 + het_het * 2 * (3/4) + het_neg

def codon_translation(rna_string):
    amino_acid_string = ""
    codons = [rna_string[i:i+3] for i in range(0,len(rna_string),3)]
    for codon in codons[:-1]:
        amino_acid_string += codon_map[codon]
    return amino_acid_string

def codon_combinations(rna_string):
    reverse_codon_count = dict()
    for acid in codon_map.values():
        if acid in reverse_codon_count:
            reverse_codon_count[acid] += 1
        else:
            reverse_codon_count[acid] = 1
    total = 1
    for char in rna_string:
        total *= reverse_codon_count[char]
    total *= reverse_codon_count['Stop']
    return total % 1_000_000

def protein_mass(protein_string):
    mass = 0
    for char in protein_string:
        mass += amino_acid_weight_map[char]
    return mass

def reverse_comp_string(dna_string):
        rev = dna_string[::-1]
        molecule_map = {"A":"T", "T":"A","C":"G","G":"C"}
        rev_comp = ""
        for char in rev:
            rev_comp += molecule_map[char]
        return rev_comp