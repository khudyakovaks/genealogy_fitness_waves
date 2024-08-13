import numpy as np
import random


def mutation(mu, population):
    """give each individual a number of new mutations distributed as Poisson(mu)"""
    N = len(population)
    mutations = np.random.poisson(mu, N)
    pop = population + mutations
    return pop


def selection_rand_drift_asexual(pop, s, ancestry, g):
    """sample next generation from the previous one and record the indices of the parents"""
    N = len(pop)
    fitnesses = mean_fitness(s, pop)  # computes fitnesses of all individuals, population mean
    # fitness and fitness variance
    probabilities = list(fitnesses / sum(fitnesses))  # probability to choose each of N individuals as parents
    parents_indices = np.random.choice(N, N, p=probabilities)  # choose N parents with replacement with probabilities
    # proportional to fitness
    next_gen = pop[parents_indices]  # children inherit parent's number of mutations
    ancestry[:, g] = parents_indices  # keep track of the identifier of the parents
    return next_gen, ancestry


def mean_fitness(s, pop):
    """"compute fitnesses of all individuals"""
    N = len(pop)
    fitnesses = (1 + s) ** (pop)  # fitness is (1+s) to the power of the number of mutations

    return fitnesses


def reconstruct_lineages(ancestry):  # reconstructs the ancestral lineages
    """"convert the list with the indices of the parents into ancestral lineages"""
    N = len(ancestry)
    n_gen = len(ancestry[0])
    lineages = np.zeros((N, n_gen))
    lineages[:, -1] = ancestry[:, -1]  # initialization
    for gen in range(n_gen - 1):  # go backward in time
        lineages[:, n_gen - gen - 2] = ancestry[
            tuple(lineages[:, n_gen - gen - 1].astype(int)), n_gen - gen - 2]  # orders lineages according to ancestry

    return lineages


def find_T2(lineages, proportion_sampled, n_gen):
    """finds average time when the MRCA of proportion_sampled*N randomly sampled pairs lived"""
    T2 = 0
    N = len(lineages)
    for i in range(int(N * proportion_sampled)):  # sample size = proportion_sampled*N
        inds = random.sample(range(N), 2)  # take two random individuals
        coinc = max(np.nonzero(lineages[inds[0]] == lineages[inds[1]])[0], default=0)  # the most recent generation
        # in which the parents of the two individuals coincide; if they don't share an ancestor, put 0

        T2 += n_gen - coinc  # count time backwards

    return T2 / int(N * proportion_sampled)  # average T2 over the sample size


def find_T3(lineages, proportion_sampled, n_gen):
    """finds average time when the MRCA of proportion_sampled*N randomly sampled triplets lived"""
    T3 = 0
    N = len(lineages)

    for i in range(int(N * proportion_sampled)):  #  sample size = N*proportion_sampled
        inds = random.sample(range(N), 3)  # sample 3 random individuals
        coinc1 = n_gen - max(np.nonzero(lineages[inds[0]] == lineages[inds[1]])[0], default=0)  # check when does
        # each possible pair has a most recent common ancestor
        coinc2 = n_gen - max(np.nonzero(lineages[inds[1]] == lineages[inds[2]])[0], default=0)
        coinc3 = n_gen - max(np.nonzero(lineages[inds[0]] == lineages[inds[2]])[0], default=0)
        T3 += max(coinc1, coinc2, coinc3)  # choose the most distant one of those

    return T3 / int(N * proportion_sampled)  # take average over the sample


def find_T4(lineages, proportion_sampled, n_gen):
    """finds average time when the MRCA of proportion_sampled*N randomly sampled quadruplets lived"""
    T4 = 0
    N = len(lineages)

    for i in range(int(N * proportion_sampled)):  #  sample size = N*proportion_sampled
        inds = random.sample(range(N), 4)  # sample 4 random individuals
        coinc1 = n_gen - max(np.nonzero(lineages[inds[0]] == lineages[inds[1]])[0], default=0)# check when does
        # each possible pair has a most recent common ancestor
        coinc2 = n_gen - max(np.nonzero(lineages[inds[0]] == lineages[inds[2]])[0], default=0)
        coinc3 = n_gen - max(np.nonzero(lineages[inds[0]] == lineages[inds[3]])[0], default=0)
        coinc4 = n_gen - max(np.nonzero(lineages[inds[1]] == lineages[inds[2]])[0], default=0)
        coinc5 = n_gen - max(np.nonzero(lineages[inds[1]] == lineages[inds[3]])[0], default=0)
        coinc6 = n_gen - max(np.nonzero(lineages[inds[2]] == lineages[inds[3]])[0], default=0)
        T4 += max(coinc1, coinc2, coinc3, coinc4, coinc5, coinc6)  # choose the most distant one of those

    return T4 / int(N * proportion_sampled)  # take average over the sample
