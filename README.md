# Coarse_Graining
Coarse graining of biochemical systems described by discrete stochastic dynamics

based on David Seiferth, Peter Sollich, and Stefan Klumpp. Phys. Rev. E 102, 062149

URL: https://link.aps.org/doi/10.1103/PhysRevE.102.062149

DOI: 10.1103/PhysRevE.102.062149

and the master thesis of the same title (see Master_Thesis.pdf).


### 'Example_Notebook.ipynb'

A few examples are given for calculating the steady-state probability distribution for master equations with finite states analytically (the analytical solution of the steady state will be compared to one obtained by a Gillespie simulation). Furthermore, this notebook shows how to calculate the entropy production (analytical mean and simulated distribution). For the molecular motor (parametrisation by Liepelt and Lipowsky 2007), the steady-state probability and entropy production are calculated. Distributions for the velocity and the entropy production of this systems are calculated from simulatios.


### 'Fig4_Distribution_Entropy_Prod_Velocity.ipynb'

Simulation results for the velocity and the entropy production of a kinesin motor. The rate constants for the six-state model are from Liepelt and Lipowsky (2007) for chemical concentrations [AT P] = [ADP] = [P] = 1 μM, stepping size l = 8 nm and an external load force F = 1 pN. 10 000 trajectories have been sampled for each model with a simulation time τ = 1200 s.
Fig. 4 refers to the PRE paper.

### 'Fig5_Iterative_Coarse_Graining/'

We iteratively eliminate states from the kinesin model to simplify the description of the molecular motor and to obtain a hierarchy of models with different levels of coarse graining. In each coarse-graining iteration, the two states that are merged are chosen such that the transition or the cycle with minimal entropy production is removed and thus the difference in entropy production between the models is minimal. In general, there will be coarse-graining steps that preserve the cycle topology and steps that change it. For the kinesin model, we can perform four iterations of the coarse-graining step. In each iteration, the model loses one state, such that we end up with a two-state network after step 4.


### 'Fig6_Fig8_Mean_and_Variances_different_models_as_function_of_force.ipynb'

To investigate whether the coarse-grained model approximates the original model well also in the presence of a load force, we plot the differences in the velocity and the entropy production between various coarse-grained models and the original model as functions of the force


### 'Dwell_time_distributions_Kinesin.ipynb'
##### refers to section B.4 in master thesis (appendix)
So far, we investigated how coarse graining affects the distribution of the velocity and the entropy production of a kinesin motor. We differed between coarse-graining mappings that preserved the network topology of the original model and between mappings that reduced the number of fundamental cycles. 
There are further quantities of interest that characterise the mechanical steps of the molecular motor. The walk of kinesin is governed by four dwell time distributions corresponding to the four possible pairs of subsequent forward and backward steps (Valleriani et al. 2008) that are called co-steps. Their distributions can be calculated from the master equation if two absorbing states are added. The effective step dynamic is non-Markovian and based on conditional mechanical steps or co-steps (Valleriani et al. 2008). 
In this section, the dwell time distributions and the steady-state probabilities for the mechanical steps and co-steps in the kinesin model of Liepelt and Lipowsky is discussed.
