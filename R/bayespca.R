#' bayespca: Regularized Principal Component Analysis via Variational Bayes inference.   
#'
#' @description A package for estimating PCA with Variational Bayes inference.   
#'
#' @details Bayesian estimation of weight vectors in PCA.
#'  To achieve regularization, the method allows specifying fixed variances
#'  in the prior distributions of the weights; alternatively, it is possible
#'  to implement Jeffrey's and  Inverse Gamma priors on such parameters.
#'  In turn, the Inverse Gamma's can have fixed shape hyperparameter; and
#'  fixed or random scale hyperparameter. Last, the method allows performing
#'  component-specific stochastic variable selection (`spike-and-slab` prior).
#'
#' 
#'
# '
#' @section Functions: 
#' \itemize{
#'     \item \code{\link{vbpca}} for model estimation; 
#'     \item \code{\link{vbpca_control}} for settings of control parameters; 
#'     \item \code{is.vbpca} for testing the class; 
#'     \item \code{\link{plothpdi}} for plotting high probability density intervals. 
#' }
#' 
#'
#' @references
#' \itemize{
#' 
#' \item [1] C. M. Bishop. 'Variational PCA'. In Proc. Ninth Int. Conf. on Artificial Neural Networks.
#' ICANN, 1999.
#'
#' \item [2] E. I. George, R. E. McCulloch (1993). 'Variable Selection via Gibbs Sampling'. 
#' Journal of the American Statistical Association (88), 881-889.
#' 
#'
#' }
#'
#' 
#' @author D. Vidotto <d.vidotto@@uvt.nl>
#'
#' @docType package  
#' @name bayespca
#' @aliases bayespca
#'
#'
#' @useDynLib bayespca, .registration = TRUE  
#' 
#' @importFrom Rcpp evalCpp
#' @importFrom ggplot2  ggplot aes geom_pointrange geom_hline coord_flip xlab 
#' @importFrom grDevices recordPlot colors 
#' @importFrom graphics plot par 
#' @importFrom stats na.omit qnorm 
NULL
