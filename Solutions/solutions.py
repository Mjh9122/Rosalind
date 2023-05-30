import numpy as np

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


def immortal_rabbit_pop_calc(months, rabbits_per_litter):
    pop, pop_past = 1, 1
    for _ in np.arange(months - 2):
        pop, pop_past = pop_past * rabbits_per_litter + pop, pop
    return pop

def mortal_rabbit_pop_calc(months, life_span):
    pass


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
