# Cox-Regression-on-the-Plane
This repository contains the R code for the simulations from the paper "Cox Regression on the Plane". The estimation is based on a generalization of the bivariate pseudo-observations approach, and includes a two-step estimation procedure for the generalized Lehmann model.

## Description
This repository contains two folders: one folder for the helper functions and another folder for the simulation functions.

### helper functions
This folder contains various helper functions tat are needed for the simulations to run.

### simulations
This folder contains three files: (1) functions to generate the simulated bivariate right censored data (corresponding to the 4 DGMs in the paper), (2) an example of the estimation procedure for the Frank NOD setting (Section 5 of our paper), and (3) a run_sim function that runs the simulations using parallel computations. Additional simulation settings can be derived similarly.
