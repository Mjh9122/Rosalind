import networkx as nx
import numpy as np
import math
import itertools

from matplotlib import pyplot as plt
from tqdm import tqdm

def immortal_rabbit_pop(months, rabbits_per_litter):
    """Population of undying rabbits after each month assuming a set rabbit per litter value

    Args:
        months (int): number of months to calculate
        rabbits_per_litter (int): number of rabbits each rabbit produces every month

    Returns:
        int: rabbit population after months months
    """
    pop, pop_past = 1, 1
    for _ in np.arange(months - 2):
        pop, pop_past = pop_past * rabbits_per_litter + pop, pop
    return pop

def mortal_rabbit_pop(months, life_span):
    """Population of rabbits assume a set life span and one rabbit per rabbit per month is born. And it takes a month to reach age of reproduction

    Args:
        months (int): number of months to calculate
        life_span (int): how many months each rabbit lives

    Returns:
        int: rabbit population are months months
    """
    pop = [0 for _ in range(life_span)]
    pop[-1] = 1
    for month in range(2, months+1):
        pop.append(sum(pop[:-1]))
        pop = pop[1:]
    return sum(pop)

def dominant_allele_by_population(homo_pos, het, homo_neg):
    """Calculate the probability of a dominant allele showing child from two randomly selected parents of the parent population given

    Args:
        homo_pos (int): homozygous positive parents
        het (int): heterozygous parents
        homo_neg (int): homozygous negative parents

    Returns:
        float: probability(not percentage) that a child of this parent population with show the dominant allele
    """
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
    """Given the number of parent matchings in a population find the expected number of dominant phenotype children

    Args:
        pos_pos (int): number of homozygous positive parent couples
        pos_het (int): number of homozygous positive - heterozygous parent couples
        pos_neg (int): number of homozygous positive - homozygous negative couples
        het_het (int): number of heterozygous positive parent couples
        het_neg (int): number of heterozygous - homozygous negative couples
        neg_neg (int): number of homozygous negative - homozygous negative couples

    Returns:
        float: expected number of positive phenotype offspring, assuming each couple has exactly two children
    """
    return pos_pos * 2 + pos_het * 2 + pos_neg * 2 + het_het * 2 * (3/4) + het_neg





