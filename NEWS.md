<<<<<<< HEAD
# bayespca 0.3.0 
* previous priors (Jeffrey's, hyperparameters Gamma prior) and SVS functionality have been restored 

=======
>>>>>>> 50009e97c685ef8e94bbfdb6fc3a466f64df3285
# bayespca 0.2.0 
* the package has been heavily revised and modified
* the SVS (Stochastic Variable Selection) functionality has been removed  (it might be re-implemented in future versions) 
* only fixed and Gamma priors for the precision parameters can be used 

# bayespca 0.1.0
* removed Type-II-Maximum-Likelihood estimation for SVS inclusion probabilities, as it was causing the algorithm to be too unstable 
* updated vignettes 
* added variable names to the model output 


# bayespca 0.0.2
* improved various computational aspects of the VB algorithm
    * more efficient SVD decomposition for large matrices 
	* bugs related to the computation of XTX fixed 
	* added print-screen iteration info 
	* fixed bug related to computation of determinants of too large posterior variance matrices (```sigmaWd``` in expectedvalues.cpp)
* added warning message for positive ```elbo``` values in case of unscaled data 
* added ```suppressWarnings``` argument to hide function warnings 
* updatetd line endings of source files to ```LF```
* added ```BugReports``` filed in description 
	


