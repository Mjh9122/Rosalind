import numpy as np


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





