#' @title Auxiliary for Controlling \code{bayespca} estimation
#'
#'
#' @aliases vbpca_control
#'
#'
#' @usage
#' vbpca_control(nstart = 1, center = TRUE, scalecorrection = 0, svdStart = TRUE,
#'      normalise = FALSE, seed = -1, plot.lowerbound = TRUE, hpdi = FALSE, probHPDI = 0.9,
#'      scaleprior = FALSE, alphatau = 0.5, betatau = 0.5, gammatau = -1, deltatau = 1,
#'      hypertype = 'common', beta1pi = 0.5, beta2pi = 0.5, v0 = 1e-04)
#'
#' @description Auxiliary control parameters for estimation of \code{vbpca} models. Internally used by
#'      \code{\link{vbpca}}, the function can be used for specifying specific parameter values.
#'
#'
#' @details
#' Use \code{seed} > 0 to let \code{\link{vbpca}} internally use the specified seed. If \code{seed} < 0, the
#' seed set in the global environment will be utilized. \cr
#' \code{scaleprior} can be set to true when the prior variance \code{tau} should be scaled by the residual variance
#' \code{sigma2}. \cr
#' \code{alphatau}, \code{betatau}, \code{gammatau}, and \code{deltatau} can be either scalar (float) values, or
#' \eqn{D} dimensional arrays, where \eqn{D} is the number of components specified in \code{\link{vbpca}}. In the
#' first case, all components will share the same hyperparameters. In the second case, each hyperparameter will be
#' component-specific. For the hyperparameters \code{gammatau} and \code{deltatau}, this type of setting can also be
#' selected with the \code{hypertype} argument. Note that the (Gamma) priors on \code{tau} will be activated only when
#' all elements of \code{gammatau} are larger than 0. \cr
#' Similarly, \code{beta1pi} and \code{beta2pi} specify the Beta prior hyperparameters for the \code{priorInclusion}
#' argument when \code{SVS = TRUE}. Both \code{beta1pi} and \code{beta2pi} can be scalar or \eqn{D} dimensional arrays.
#' The Beta hyperpriors will be activated only when all elements of \code{beta1pi} are
#' larger than 0.
#'
#'
#' @param nstart integer; \cr
#'               number of (sequential) estimation of the \code{vbpca} algorithm.
#'
#' @param center bool; \cr
#'               boolean indicating whether to center the variables of the input data before model estimation.
#'
#' @param scalecorrection integer; \cr
#'                        factor used for the scaling of the variables prior to the estimation step.
#'                        When smaller than 0, no scaling will be performed.
#'
#' @param svdStart bool; \cr
#'                 boolean indicating whether the values of the SVD decomposition of the input matrix shoul be
#'                 used as starting values.
#'
#' @param normalise bool; \cr
#'                  logical argument indicating whether the elements of the weight matrix \eqn{W} should be
#'                  normalised (after the estimation of the model).
#'
#' @param seed integer; \cr
#'             seed used for the random initialization of the model.
#'
#' @param plot.lowerbound bool; \cr
#'                        boolean indicating whether the function should plot the traceplot of the ELBO values
#'                        calculated across the Variational Bayes iterations.
#'
#' @param hpdi bool; \cr
#'             boolean denoting whether the high posterior density intervals of \eqn{W} should be computed.
#'
#' @param probHPDI float; \cr
#'                 the desired probability level of the HPD intervals.
#'
#' @param scaleprior bool; \cr
#'                   logical; when TRUE, the prior variances \eqn{tau} will be scaled by the residual variance \code{sigma2}.
#'
#' @param alphatau float or array_like; \cr
#'                 shape parameter for the Inverse Gamma hyperpriors. It can be scalar, or \eqn{D} dimensional.
#'
#' @param betatau float or array_like; \cr
#'                scale parameter for the Inverse Gamma hyperprior. It can be scalar, or \eqn{D} dimensional.
#'
#' @param gammatau float or array_like; \cr
#'                 shape parameter for the Gamma hyperprior on \code{betatau}. It can be scalar, or \eqn{D} dimensional.
#'                 In order for the Gamma hyperprior to be activated, all elements of \code{gammatau} must be larger than 0.
#'
#' @param deltatau float or array_like; \cr
#'                 scale parameter for the Gamma hyperprior on \code{betatau}. It can be scalar, or \eqn{D} dimensional.
#'
#' @param hypertype character; \cr
#'                  string denoting the type of update for the hyperparameters of the Gamma prior on \code{betatau}.
#'                  Use \code{'global'} for sharing information among all components; use \code{'component'} for
#'                  component-specific updates; use \code{'local'} for element-specific updates.
#'
#' @param beta1pi float or array_like; \cr
#'                shape parameter which models the probability of inclusion in the Beta hyperprior (\code{SVS = TRUE}).
#'                It can be a scalar, or a \eqn{D} dimensional vector. When at least one element of \code{beta1pi}
#'                is smaller than 0, no updates will be performed on \code{priorInclusion}. When all elements are
#'                equal to 0, or all elements are positive and at least one is equal to 0, type-II ML updates will occur.
#'                When all the elements are larger than 0, a Beta hyperprior is assumed on \code{priorInclusion}.
#'
#' @param beta2pi float or array_like; \cr
#'                shape parameter which models the probability of exclusion in the Beta hyperprior (\code{SVS = TRUE}).
#'                It can be a scalar, or a \eqn{D} dimensional vector.
#'
#' @param v0 float; \cr
#'           scalar value between 0 and 1 (possibly close to 0) which specifies how much the prior variance of the `slab`
#'           variance \code{tau} should be scaled for the `spike` component of the prior.
#'
#'
#' @return a list containing the control parameters specified by the user, as well as the unspecified default values.
#'
#'
#' @author D. Vidotto <d.vidotto@@uvt.nl>
#'
#' @seealso \code{\link{vbpca}}
#'
#' @examples
#'
#' # Specify controls for Inverse Gamma(1, .01) prior for W;
#' # and Beta(5, 1) prior for priorInclusion
#' ctrl <- vbpca_control(alphatau = 1, betatau = .01, beta1pi = 5, beta2pi = 1)
#'
#'
#'
#'



#' @export
vbpca_control <- function(nstart = 1, center = TRUE, scalecorrection = 0, svdStart = TRUE, normalise = FALSE, seed = -1, plot.lowerbound = TRUE, hpdi = FALSE, probHPDI = 0.9, 
    scaleprior = FALSE, alphatau = 0.5, betatau = 0.5, gammatau = -1, deltatau = 1, hypertype = "common", beta1pi = 0.5, beta2pi = 0.5, v0 = 1e-04) {
    
    controlList <- list(nstart = nstart, center = center, scalecorrection = scalecorrection, svdStart = svdStart, normalise = normalise, seed = seed, plot.lowerbound = plot.lowerbound, 
        hpdi = hpdi, probHPDI = probHPDI, scaleprior = scaleprior, alphatau = alphatau, betatau = betatau, gammatau = gammatau, deltatau = deltatau, hypertype = hypertype, 
        beta1pi = beta1pi, beta2pi = beta2pi, v0 = v0)
    
    return(controlList)
    
}
