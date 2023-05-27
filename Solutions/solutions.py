import numpy as np

def rabbit_pop_calc(months, rabbits_per_litter):
    pop, pop_past = 1, 1
    for _ in np.arange(months-2):
        pop, pop_past = pop_past*rabbits_per_litter+pop, pop
    return pop
