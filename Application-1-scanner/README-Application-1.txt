### Instructions for the covariance estimation method for Application 1 ###

Before using the covariance estimation method, we recommend exploring code that was provided in a previous publication that forms the basis for the new approach in this publication: 
- Estimating model parameters via Particle Metropolis within Gibbs (PMwG) algorithm (osf.io/5b4w3).

* LBA_PMwG_v1.m : This is the primary file to run the method. Start here. This file loads the data as required:
   -> LBA_Forstmann_scanner.mat : Data from Forstmann et al. (2008) in format expected by the algorithm, for sessions in and out of the scanner. 

When applying the covariance estimation method to different data sets, the following scripts will require editing depending on the design of the new data set and the model that is estimated. The main file (above) will also require editing. 
* LBA_MC_v1.m
* LBA_CMC_v1.m
* reshape_A.m
* reshape_b.m
* reshape_v.m
* reshape_tau.m
* plot_result.m
* compute_correlation.m

The remaining files are subsidiary Matlab functions. These don't require editing for application to different data sets.
