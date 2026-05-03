# Cox-Regression-on-the-Plane

This repository contains the R code for the simulation study in the paper
**“Cox Regression on the Plane.”**
A preprint is available at: https://arxiv.org/abs/2509.12473

The implementation is based on a generalization of the bivariate pseudo-observations approach and includes a two-step estimation procedure for the generalized Lehmann model.

---

## 📦 Overview

This repository allows you to:

* Simulate bivariate survival data under multiple data-generating mechanisms (DGMs)
* Fit the simple Lehmann model and the generalized Lehmann model
* Apply the proposed two-step estimation procedure
* Compare with copula-based models
* Reproduce representative simulation results reported in the paper

The generalized Lehmann model is estimated via a two-step procedure: marginal parameters are estimated first, followed by estimation of the dependence parameter conditional on the first step. The example script provides a minimal working pipeline from data generation to estimation.

---

## Repository Structure

The repository contains two main folders:

### `helper_functions`

This folder contains auxiliary functions required for the simulations, including routines for estimation, variance calculation, and supporting computations.

### `simulations`

This folder contains:

1. **Data generation functions**
   Functions to simulate bivariate right-censored data under the four DGMs considered in the paper.

2. **Example script**
   An example illustrating the estimation procedure for the Frank NOD setting (see Section 5 of the paper).

3. **Simulation runner**
   A `run_sim` function that executes simulation studies using parallel computation.

Additional simulation scenarios can be implemented by modifying the example script accordingly.

---

## Reproducibility Notes

* The scripts are designed to reproduce the main components of the simulation study.
* Simulation settings (e.g., sample size, number of replications, DGMs) can be adjusted within the scripts.
* Runtime depends on the chosen settings; larger simulation runs may be computationally intensive and may benefit from parallel or server-based execution.

---

## Citation

If you use this code, please cite the corresponding paper:

*Travis-Lumer, Y., Mandel, M., Fabian I. D., Betensky, R. A., and Gorfine, M. (2025).*
*Cox Regression on the Plane.*
