# Gibbs sampler for dynamic community detection

This repository contains the source code for a simulation study in my master's thesis 
Consistent Bayesian community detection (Aalto University, 2021).
The simulation study examines the effect of the number of observed
network layers on the classification accuracy of a community recovery algorithm. 
The Gibbs sampler presented here is an extension of
that of Nowicki and Snijders ([2001](#sources)), which takes
as an input a tensor of independent and identically distributed adjacency matrices.

The sampler's implementation follows closely that of MFM-SBM of Geng et al. ([2019](#sources)), depicted fully in 
the associated [supplemental material](https://www.tandfonline.com/doi/suppl/10.1080/01621459.2018.1458618).

Built on R version 4.1.0 (2021-05-18).

## Sources: {#sources}
Nowicki, K. and Snijders, T. A. B. (2001), ‘Estimation and prediction
for stochastic blockstructures’, Journal of the American Statistical Association
96(455), 1077–1087.

Geng, J., Bhattacharya, A. and Pati, D. (2019), ‘Probabilistic community detection with unknown number of communities’, Journal of the American Statistical
Association 114(526), 893–905.

