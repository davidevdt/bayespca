# bayespca
A package for regularized Principal Component Analysis via Variational Bayes methods. 

## Author
Davide Vidotto <d.vidotto@uvt.nl> 

## Description
`bayespca` performs Bayesian estimation of weight vectors in PCA.
    To achieve regularization, the method allows specifying fixed variances
    in the prior distributions of the weights; alternatively, it is possible
    to implement Jeffrey's and  Inverse Gamma priors on such parameters.
    In turn, the Inverse Gamma's can have fixed shape hyperparameter; and
    fixed or random scale hyperparameter. Last, the method allows performing
    component-specific Stochastic Variable Selection (`spike-and-slab` prior).
    
    Check the ```vignettes``` and package documentation for further details. 

## Functions

* ```vbpca``` for model estimation;
* ```vbpca_control``` for settings of control parameters; 
* ```is.vbpca``` for testing the class; 
* ```plothpdi``` for plotting high probability density intervals. 

## Version
0.0.1

## Depends 
R (>= 3.3.3)

## License 
GPL-2



