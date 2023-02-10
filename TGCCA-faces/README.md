# TGCCA-faces
R code to evaluate TGCCA on real data using the Multi-PIE Face dataset [1].

The data we used has been preprocessed for the work of Tian et al. (2018) [2] and is
available on their [github repository](https://github.com/bluer555/CR-GAN/blob/master/README.md).
We take the first 100 subjects to form the training set and the next 100 to form the testing set.
The goal of this experiment is to learn a latent subspace between two views of the training
set and use this subspace to pair images of subjects from the testing set in these two views.
15 illumination conditions are used in the training set. 1 to 15 illumination conditions
are used in the testing set.

## Prerequisites
- R package RGCCA. Instructions to install it are in the main folder.
- R packages dplyr, tidyr, imager, abind, ggplot2, viridis and lpSolve.

## Run and compare models
Models are run only once but the sampling of the illumination conditions
is repeated 100 times. The constructed graph reports the median matching
accuracies and the intervals containing 95% of the achieved accuracies.
The figure reported in the paper was produced running the following command: \
`Rscript path_to_TGCCA-faces/evaluation.R path_to_TGCCA-faces`


## References
1. Gross, R., Matthews, I., Cohn, J., Kanade, T., and Baker, S. (2008). Multi-pie. In 2008 8th IEEE International Conference on Automatic Face & Gesture Recognition, pages 1–8.
2. Tian, Y., Peng, X., Zhao, L., Zhang, S., and Metaxas, D. N. (2018). Cr-gan: Learning complete representations for multi-view generation.
