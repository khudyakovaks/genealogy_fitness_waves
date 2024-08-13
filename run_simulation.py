from functions import *
from parameters import *
import numpy as np


T2s = np.zeros((len(ses), len(mus)))
T3s = np.zeros((len(ses), len(mus)))
T4s = np.zeros((len(ses), len(mus)))

for i in range(len(ses)):  # for each s
    for j in range(len(mus)):  # for each mu
        print("s ", i, " out of ", len(ses), "; mu ", j, " out of ", len(mus))
        s = ses[i]
        mu = mus[j]
        for r in range(n_replicates):  # for each replicate
            population = np.zeros(N)  # initialize a population
            ancestry = np.zeros((N, n_gen))  # initialize a list of parent indices for all generations
            for g in range(n_gen):  # for each generation
                pop_mut = mutation(mu, population)  # population after mutation
                population, ancestry = selection_rand_drift_asexual(pop_mut, s, ancestry, g)  # population after
                # selection
            lineages = reconstruct_lineages(ancestry)  # reconstructs ancestral lineages from list of parent indices
            T2s[i, j] += find_T2(lineages, proportion_sampled, n_gen) / n_replicates  # record T2 averaged across
            # replicates
            T3s[i, j] += find_T3(lineages, proportion_sampled, n_gen) / n_replicates  # record T3 averaged across
            # replicates
            T4s[i, j] += find_T4(lineages, proportion_sampled, n_gen) / n_replicates  # record T4 averaged across
            # replicates

np.save("T2s_N{}_n_gen{}_n_replicates{}_proportion_sampled_{}".format(N, n_gen, n_replicates, proportion_sampled), T2s)
np.save("T3s_N{}_n_gen{}_n_replicates{}_proportion_sampled_{}".format(N, n_gen, n_replicates, proportion_sampled), T3s)
np.save("T4s_N{}_n_gen{}_n_replicates{}_proportion_sampled_{}".format(N, n_gen, n_replicates, proportion_sampled), T4s)
