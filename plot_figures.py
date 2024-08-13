from parameters import *
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


plt.rcParams.update({'font.size': 22})
colormap = sns.color_palette("coolwarm", n_colors=10)


# 1. load the output of 'run_simulation.py'
T2s = np.load(
    "T2s_N{}_n_gen{}_n_replicates{}_proportion_sampled_{}.npy".format(N, n_gen, n_replicates, proportion_sampled))
T3s = np.load(
    "T3s_N{}_n_gen{}_n_replicates{}_proportion_sampled_{}.npy".format(N, n_gen, n_replicates, proportion_sampled))
T4s = np.load(
    "T4s_N{}_n_gen{}_n_replicates{}_proportion_sampled_{}.npy".format(N, n_gen, n_replicates, proportion_sampled))

# 2. compute the scaling of T2 in the background selection regime
ces = mus / (-ses[:, np.newaxis] * np.log(
    N))  # constant for the scaling from Igelbrink et al. 2023, as explained in the main text
ces[ces > 1] = 1
ces[-1] = np.array([0] * len(mus))  # the scaling factor in the case when s=0

T2s_scaled = T2s / N ** (1 - ces)
T2s_scaled[T2s_scaled > 2] = 2  # cut values above 2, as explained in the main text
T2_scaled_df = pd.DataFrame(T2s_scaled, columns=np.around(mus, 3), index=np.around(ses, 3))

# 4. compute log T2 / log N
logT2_over_logN = pd.DataFrame(np.log(T2s) / np.log(N), columns=np.around(mus, 3), index=np.around(ses, 3))

# 5. infer alpha (parameter of the Beta coalescent)
alphas = (2 * (1 - (T3s / T2s))) / (2 * (T3s / T2s) - 3)  # values of alpha inferred from <T3>/<T2>
alphas[alphas > 2] = 2
alphas_dataframe = pd.DataFrame(alphas, columns=np.around(mus, 3), index=np.around(ses, 3))

x = T4s / T2s
alphas_T4T2 = (14 - 9 * x - np.sqrt(76 - 60 * x + 9 * x ** 2)) / (2 * (-5 + 3 * x))  # same, but inferred from <T4>/<T2>
alphas_T4T2[alphas_T4T2 > 2] = 2
alphas_T4T2_dataframe = pd.DataFrame(alphas_T4T2, columns=np.around(mus, 3), index=np.around(ses, 3))

# 6. get values of s that generate Г = 1.2, Г = 1 as a function of mu to plot Г-lines
ses_to_give_gammas_one_half = 2 * mus / (np.log(mus * N))  # Г = 1/2 line
ses_to_give_gammas_one = mus / (np.log(mus * N))  # Г = 1 line

# 7. plot the heatmap of alpha values inferred from <T3>/<T2>, together with two critical Г lines
fig, ax = plt.subplots(figsize=(20, 12))
ax = sns.heatmap(alphas_dataframe, linewidth=0.5, cmap=colormap, center=1.1, robust=True)
ax2 = plt.twinx()
sns.lineplot(data=-ses_to_give_gammas_one_half, linewidth=5, color='black', linestyle='dashed', ax=ax2)
sns.lineplot(data=-ses_to_give_gammas_one, linewidth=5, color='black', linestyle='dashed', ax=ax2)
plt.ylim(ses[-1], ses[0])
ax2.set_yticks([])
ax.set_ylabel(r'$s$', fontsize=35)
ax.set_xlabel(r'$\mu$', fontsize=35)
plt.title(r"Values of $\alpha$ from $<T_3> / <T_2>, N = %i$" % N, fontsize=35)
plt.savefig(
    'Heatmap_T3T2_to_alphas_mu_vs_s_N{}_n_gen{}_n_replicates{}_proportion_sampled_{}.png'.format(N, n_gen, n_replicates,
                                                                                                 proportion_sampled))

# plot the heatmap of alpha values inferred from <T4>/<T2>, together with two critical Г lines
fig, ax = plt.subplots(figsize=(20, 12))
ax = sns.heatmap(alphas_T4T2_dataframe, linewidth=0.5, cmap=colormap, center=1.1, robust=True)
ax2 = plt.twinx()
sns.lineplot(data=-ses_to_give_gammas_one_half, linewidth=5, color='black', linestyle='dashed', ax=ax2)
sns.lineplot(data=-ses_to_give_gammas_one, linewidth=5, color='black', linestyle='dashed', ax=ax2)
plt.ylim(ses[-1], ses[0])
ax2.set_yticks([])
ax.set_ylabel(r'$s$', fontsize=35)
ax.set_xlabel(r'$\mu$', fontsize=35)
plt.title(r"Values of $\alpha$ from $<T_4> / <T_2>, N = %i$" % N, fontsize=35)
plt.savefig(
    'Heatmap_T4T2_to_alphas_both_gamma_mu_vs_s_N{}_n_gen{}_n_replicates{}_proportion_sampled_{}.png'.format(N, n_gen,
                                                                                                            n_replicates,
                                                                                                            proportion_sampled))

# plot the heatmap of <T2>/N^{1-c}
fig, ax = plt.subplots(figsize=(20, 12))
ax = sns.heatmap(T2_scaled_df, linewidth=0.5, cmap=colormap)
ax2 = plt.twinx()
sns.lineplot(data=-ses_to_give_gammas_one_half, linewidth=5, color='black', linestyle='dashed', ax=ax2)
sns.lineplot(data=-ses_to_give_gammas_one, linewidth=5, color='black', linestyle='dashed', ax=ax2)
plt.ylim(ses[-1], ses[0])
ax2.set_yticks([])
ax.set_ylabel(r'$s$', fontsize=35)
ax.set_xlabel(r'$\mu$', fontsize=35)
plt.title(r"$\frac{<T_2>}{N^{1-c}}, N = %i$" % N, fontsize=35)
plt.savefig(
    'Heatmap_scaled_T2_mu_vs_s_N{}_n_gen{}_n_replicates{}_proportion_sampled_{}.png'.format(N, n_gen, n_replicates,
                                                                                            proportion_sampled))

# plot the heatmap of log <T2> / log N
fig, ax = plt.subplots(figsize=(20, 12))
ax = sns.heatmap(logT2_over_logN, linewidth=0.5, cmap=colormap)
ax2 = plt.twinx()
sns.lineplot(data=-ses_to_give_gammas_one_half, linewidth=5, color='black', linestyle='dashed', ax=ax2)
sns.lineplot(data=-ses_to_give_gammas_one, linewidth=5, color='black', linestyle='dashed', ax=ax2)
plt.ylim(ses[-1], ses[0])
ax2.set_yticks([])
ax.set_ylabel(r'$s$', fontsize=35)
ax.set_xlabel(r'$\mu$', fontsize=35)
plt.title(r"$\log(<T_2>) / \log(N), N = %i$" % N, fontsize=35)
plt.savefig('Heatmap_logT2_over_logN_mu_vs_s_N{}_n_gen{}_n_replicates{}_proportion_sampled_{}.png'.format(N, n_gen,
                                                                                                          n_replicates,
                                                                                                          proportion_sampled))
