# TGCCA-2D_simulations
R code to generate simulations to test TGCCA and compare it with TCCA  [1] and [2DCCA](https://github.com/youlinchen/TCCA) [2].

## Prerequisites
- R package RGCCA. Instructions to install it are in the main folder.
- R packages reticulate, multiway, R.matlab and abind.

## Generate data
To generate data, the script `generate.R` is used. It generates 100 folds of 100 observations each from the specified shapes, noise structure, signal-to-noise ratio (SNR) and correlations.
The data used in the paper was generated running the following command: \
`Rscript path_to_TGCCA-2D_simulations/generate.R 100 100 SNR 0.8,0.8,0.8,0.8,0.8 square,gas,cross,cross_little,vector path_to_TGCCA-2D_simulations information,parking,restaurant,cup,null`\
for values of SNR `0.1`, `0.3`, `0.5`, and `1`.

## Run and compare models
To help running and comparing different models, the script `make_scripts.R` has been created. \
To generate the scripts that need to be run to reproduce the experiments, the following commands have been run
with `n = 100` and `N = 10000`, `n = 200` and `N = 10000`, `n = 300` and `N = 9900`, `n = 500` and `N = 10000`, and `n = 1000` and `N = 10000`;
and `experiment = 2B` or `experiment = 5B`: \
`Rscript path_to_TGCCA-2D_simulations/make_scripts.R n N path_to_TGCCA-2D_simulations experiment`

Running the script will create commands, in folder `cmd`, to run in a shell, that will conduct the analysis.\
Results will appear in folder `results` and final comparisons will appear in folder `comparisons`. A figure summarizing the comparison for the different SNR levels
as well as a table per SNR level are generated.

## Additional analysis
The additional analysis shown in the supplementary material about the impact of the initialization point and the contributions of the different rank-1 factors of
spTGCCA were respectively done using the scripts `try_initialization.R` and `compare_weights.R`.

## References
1. Chen, Y.-L., Kolar, M., and Tsay, R. S. (2021). Tensor canonical correlation analysis with convergence and statistical guarantees. Journal of Computational and Graphical Statistics, 0(0):1–17.
2. Min, E. J., Chi, E. C., and Zhou, H. (2019). Tensor canonical correlation analysis.Stat, 8(1):e253. e253 sta4.253.
