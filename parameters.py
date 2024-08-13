"""Parameters of the simulation to be configured."""

import numpy as np

N = int(100)  # population size
ses = np.arange(-0.15, 0.0025, 0.0025)  # list of selection coefficients
mus = np.arange(0.01, 0.3, 0.005)  # list of mutation rates

n_gen = int(5*N)  # number of generations
n_replicates = 100  # number of independent replicates
proportion_sampled = 1  # proportion_sampled*N = number of pairs, triplets and quadruples to compute <T2>, <T3> and <T4>
