# bayespca 0.0.2

* improved various computational aspects of the VB algorithm
    * more efficient SVD decomposition for large matrices	implemented 
	* bugs related to the computation of XTX fixed 
	* added print-screen iteration info 
	* fixed bug related to computation of determinants of too large posterior variance matrices (```sigmaWd``` in expectedvalues.cpp)
* added warning message for positive ```elbo``` values in case of unscaled data 
* added ```suppressWarnings``` argument to hide function warnings 
	