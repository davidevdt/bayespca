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
	


# bayespca 0.1.0
* removed Type-II-Maximum-Likelihood estimation for SVS inclusion probabilities, as it was causing the algorithm to be too unstable 
* updated vignettes 
* added variable names to the model output 