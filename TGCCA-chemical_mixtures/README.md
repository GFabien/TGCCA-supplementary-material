# TGCCA-chemical_mixtures
R code to evaluate TGCCA on real data and compare it with CMTF from E. Acar, T. G. Kolda, and D. M. Dunlavy and ACMTF
from E. Acar, E. E. Papalexakis, G. Gurdeniz,  M. A. Rasmussen,  A. J. Lawaetz, M. Nilsson and  R. Bro.

The data is the data used to validate ACMTF and is available from this [link](http://www.models.life.ku.dk/joda/prototype).
As in the original paper describing ACMTF, we use the NMR (Nuclear Magnetic Resonance)
spectroscopy and LC-MS (Liquid Chromatography - Mass Spectrometry) blocks.
From these two blocks, we aim to retrieve the concentrations of the 5 chemicals
used to create the 28 mixtures composing the data.

## Prerequisites
- R package RGCCA. Instructions to install it are in the main folder.
- R packages mutliway, R.matlab, ggplot2, viridis and abind.
- Matlab, Matlab [tensor toolbox](https://www.tensortoolbox.org/), Matlab [CMTF toolbox](http://www.models.life.ku.dk/~acare/CMTF_Toolbox.html), and Matlab [Poblano toolbox](https://github.com/sandialabs/poblano_toolbox).

## Run and compare models
All models were run 100 times with random initialization points. The reported
results include means and standard deviations of the cosines between the
true chemical concentrations and the estimated ones, the cosines obtained
by the models with the best criteria over the multiple runs and the computation
times (means and standard deviations).

### Run CMTF and ACMTF
CMTF and ACMTF models were run using the `run_cmtf_and_acmtf.m` file. To directly
run it from command line, you can do \
`matlab  -nodisplay -nosplash -nodesktop -r "data_path='path_to_TGCCA-chemical_mixtures';run('path_to_TGCCA-chemical_mixtures/run_cmtf_and_acmtf.m');exit()"`

### Run TGCCA
Two versions of TGCCA were run, one with a first component of rank 2 and the
other one with all components having rank 1. This was done using the
`tgcca_analysis.R` file. To run it from command line, you can do \
`Rscript path_to_TGCCA-chemical_mixtures/tgcca_analysis.R path_to_TGCCA-chemical_mixtures`

### Compile results
To analyse the results, produce table and plots, the `compile_results.R` file
was used: \
`Rscript path_to_TGCCA-chemical_mixtures/compile_results.R path_to_TGCCA-chemical_mixtures`

## References
1. E. Acar, T. G. Kolda, and D. M. Dunlavy, All-at-once Optimization for Coupled Matrix and Tensor Factorizations, KDD Workshop on Mining and Learning with Graphs, 2011 (arXiv:1105.3422v1).
2. E. Acar, E. E. Papalexakis, G. Gurdeniz,  M. A. Rasmussen,  A. J. Lawaetz, M. Nilsson,  R. Bro, Structure-Revealing Data Fusion, BMC Bioinformatics, 15: 239, 2014.
