# TGCCA-3D_simulations
R code to generate simulations to test TGCCA and compare it with TCCA [1] on simulations with 3D shapes for canonical vectors.

## Prerequisites
- R package RGCCA. Instructions to install it are in the main folder.
- R packages multiway, R.matlab and abind.

## Generate data
To generate data, the script `generate.R` is used. It generates 100 folds of 100 observations each from the specified signal-to-noise ratio (SNR) and correlations.
The data used in the paper was generated running the following command: \
`Rscript path_to_TGCCA-3D_simulations/generate.R 100 SNR 0.8,0.8,0.8,0.8,0.8 path_to_TGCCA-3D_simulations`\
for values of SNR `0.1`, `0.3`, `0.5`, and `1`.

## Run and compare models
To help running and comparing different models, the script `make_scripts.R` has been created. \
To generate the scripts that need to be run to reproduce the experiments, the following commands have been run
with `n = 100` and `N = 10000`, `n = 200` and `N = 10000`, `n = 500` and `N = 10000`, and `n = 1000` and `N = 10000`;
and `experiment = 2B` or `experiment = 5B`: \
`Rscript path_to_TGCCA-3D_simulations/make_scripts.R n N path_to_TGCCA-3D_simulations experiment`

Running the script will create commands, in folder `cmd`, to run in a shell, that will conduct the analysis.\
Results will appear in folder `results` and final comparisons will appear in folder `comparisons`. A figure summarizing the comparison for the different SNR levels
as well as a table per SNR level are generated.

## References
1. Chen, Y.-L., Kolar, M., and Tsay, R. S. (2021). Tensor canonical correlation analysis with convergence and statistical guarantees. Journal of Computational and Graphical Statistics, 0(0):1–17.
