# genealogy_fitness_waves
This repository contains the code accompanying Khudiakova, Boenkost, and Tourniaire (2024) bioRxiv.

We follow the genealogies in the individual-based Wright-Fisher simulation of an asexual population subject to deleterious mutations and multiplicative selection. 

`parameters.py` contains configurable simulation parameters.

`run_simulation.py` contains the code to run the simulations.

`plot_figures.py` produces the figures.

`functions.py` defines all the functions necessary to run the simulation.

Keeping track of the genealogies is computationally heavy, and to perform the run for large N and many replicates (as in the manuscript), the user might need to set up a parallelization on a high-performance cluster.

## Requirements and Installation:

**Requirements**

- Python 3.11 or higher
- Packages: numpy, seaborn, pandas, matplotlib

**Installation:**

`git clone https://github.com/khudyakovaks/genealogy_fitness_waves`

## Usage:
**1) Specify the parameters**

The user needs to specify the following parameters of the simulation by editing `parameters.py`: 
  * `N` - population size
  * `ses` - a list of selection coefficients
  * `mus` - a list of mutation rates
  
There is an option to specify additional parameters in the same file:
  * `n_replicates` - number of independent runs of the simulation. Increasing the number of replicates decreases the noisiness of the outcome, but increases the run time
  * `n_gen` - number of generations for the forward simulation. We recommend setting it to `5N` generations to ensure that all lineages present in the last generation originate from one ancestral individual
  * `proportion_sampled` - sample size to compute <T_k> = proportion_sampled*N. The genealogy trees are highly noisy, so sampling too many times from the same tree biases the statistics. At the same time, sampling from the same tree several times allows to reduce the number of independent runs. Through trial and error, we found that setting this parameter to `1` (i.e., sampling `N` pairs, triplets, and quadruples) leads to good results.

**2) Run `run_simulation.py`**

The output of the run is three `*.npy` files: the estimated average time to coalescence of 2, 3, and 4 individuals.

**3) Plot the figures**

After `run_simulation.py` is finished, the user can run `plot_figures.py` to produce the heatmaps as in the manuscript.
