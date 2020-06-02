#' bayespca: Regularized Principal Component Analysis via Variational Bayes inference.
#'
#' @description A package for estimating PCA with Variational Bayes inference.
#'
#' @details Bayesian estimation of weight vectors in PCA.
#'  To achieve regularization, the method allows specifying fixed precisions
#'  in the prior distributions of the weights; alternatively, it is possible
#'  to specify Gamma priors for such parameters. The method allows
#'  for variable selection through Automatic Relevance Determination.
#'
#'
# '
#' @section Functions:
#' \itemize{
#'     \item \code{\link{vbpca}} for model estimation;
#'     \item \code{is.vbpca} for testing the class;
#'     \item \code{\link{plotheatmap}} for plotting the precision and weights matrices;
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
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_hline geom_tile coord_flip xlab ylab scale_fill_gradient2 ggtitle  theme_bw theme scale_fill_distiller element_blank element_text 
#' @importFrom scales squish
#' @importFrom grDevices recordPlot colors
#' @importFrom graphics plot par
#' @importFrom stats na.omit qnorm
NULL
