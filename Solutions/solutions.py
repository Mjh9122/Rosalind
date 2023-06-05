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

def blossom_map():
    with open('Solutions/blossum62.txt') as f:
        lines = f.readlines()
        matrix = [line.strip().split() for line in lines]
        top_labels = matrix[0]
        matrix = matrix[1:]
        side_labels = [matrix[i][0] for i in range(len(matrix))]        
        matrix = [matrix[i][1:] for i in range(len(matrix))]
        blossom_map = dict()
        for i, top in enumerate(top_labels):
            for j, side in enumerate(side_labels):
                blossom_map[top+side] = int(matrix[j][i])

        return blossom_map

def pam_map():
    with open('Solutions/pam50.txt') as f:
        lines = f.readlines()
        matrix = [line.strip().split() for line in lines]
        top_labels = matrix[0]
        matrix = matrix[1:]
        side_labels = [matrix[i][0] for i in range(len(matrix))]        
        matrix = [matrix[i][1:] for i in range(len(matrix))]
        blossom_map = dict()
        for i, top in enumerate(top_labels):
            for j, side in enumerate(side_labels):
                blossom_map[top+side] = int(matrix[j][i])

        return blossom_map
    
def mismatch_map(match_score, missmatch_score):
    with open('Solutions/blossum62.txt') as f:
        lines = f.readlines()
        matrix = [line.strip().split() for line in lines]
        top_labels = matrix[0]
        side_labels = [matrix[i][0] for i in range(len(matrix))]        
        mismatch_map = dict()
        for i, top in enumerate(top_labels):
            for j, side in enumerate(side_labels):
                if top == side:
                    mismatch_map[top+side] = 1
                else:
                    mismatch_map[top+side] = -1
        return mismatch_map


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
        try:
            amino_acid_string += codon_map[codon]
        except:
            return ""
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

def to_rna(dna_string):
        return dna_string.replace('T', 'U')

def edit_distance_weighted(string1, string2, map_func, gap_pen):
        mapper = map_func()
        rows = string1
        cols = string2
        dist = [[0 for _ in range(len(cols)+1)] for _ in range(len(rows)+1)]
        for r in range(1, len(rows)+1):
            for c in range(1, len(cols)+1):
                dist[r][c] = max(dist[r-1][c] - gap_pen, dist[r-1][c-1] + mapper[rows[r-1]+cols[c-1]], dist[r][c-1] - gap_pen, 0)
        return dist[-1][-1]

def edit_distance(dna_string1, dna_string2):
        rows = dna_string1
        cols = dna_string2
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

