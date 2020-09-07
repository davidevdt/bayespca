# bayespca: Regularized Principal Component Analysis via Variational Bayes inference   
An R package for regularized Principal Component Analysis via Variational Bayes methods.

## Author
Davide Vidotto <d.vidotto@uvt.nl>

## Description
`bayespca` performs Bayesian estimation of weight vectors in PCA.
    To achieve regularization, the method allows specifying fixed precisions
    in the prior distributions of the weights; alternatively, it is possible
    to implement Gamma priors on such parameters. The method allows
    for variable selection through Automatic Relevance Determination.
    Check the ```vignettes``` and package documentation for further details.

## Functions

* ```vbpca``` for model estimation
* ```vbpca_control``` for settings of control parameters
* ```is.vbpca``` for testing the class
* ```plotheatmap``` for plotting the precision and weights matrices;
* ```plothpdi``` for plotting high probability density intervals

## Install
devtools::install_github("davidevdt/bayespca")

## Version
0.3.0

## Depends
R (>= 3.3.3)

## License
GPL-2
