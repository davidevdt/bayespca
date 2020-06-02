#' @title Regularized Variational Bayes Principal Compnent Analysis (vbpca).
#'
#'
#' @aliases vbpca.default print.vbpca summary.vbpca is.vbpca print.summary.vbpca
#'
#'
#' @usage
#' vbpca(X, D = 1, nstart = 1, maxIter = 500, tolerance = 1e-05, center = TRUE,
#' 		 scalecorrection = 1, svdStart = TRUE, verbose = FALSE, 
#' 		normalise = FALSE, seed = 1, tau = 1, updatetau = FALSE, 
#' 		alphatau = 0, betatau = 0, plot.lowerbound = TRUE, 
#' 		hpdi = FALSE, probHPDI = 0.9, global.var = FALSE, 
#' 		suppressWarnings = FALSE)
#'
#' \method{vbpca}{default}(X, D = 1, nstart = 1, maxIter = 500, tolerance = 1e-05, center = TRUE,
#' 		 scalecorrection = 1, svdStart = TRUE, verbose = FALSE, 
#' 		normalise = FALSE, seed = 1, tau = 1, updatetau = FALSE, 
#' 		alphatau = 0, betatau = 0, plot.lowerbound = TRUE, 
#' 		hpdi = FALSE, probHPDI = 0.9, global.var = FALSE, 
#' 		suppressWarnings = FALSE)
#'
#' \method{print}{vbpca}(x, ...)
#'
#' \method{summary}{vbpca}(object, ...)
#'
#' is.vbpca(object)
#'
#'
#' @description Estimation of regularized PCA with a Variational Bayes algorithm.
#'
#' @details
#' The function allows performing PCA decomposition of an \eqn{(I, J)} input matrix \eqn{X}.
#' For D principal components, the factorization occurs through:
#'
#' \deqn{ X = X W P^T + E }
#'
#' where \eqn{P} is the \eqn{(J, D)} orthogonal loading matrix (\eqn{P^T P = I}) and \eqn{W} is the
#' \eqn{(J, D)} weight matrix. E is an \eqn{(I, J)} residual matrix.
#' Principal components are defined by \eqn{X W}. In this context, focus of the inference is on the
#' weight matrix \eqn{W}. The Variational Bayes algorithm treats the elements
#' of \eqn{W} as latent variables; \eqn{P} and \code{sigma^2} (the variance of the residuals) are
#' fixed parameters instead.
#'
#' In order to regularize the elements of \eqn{W}, a Multivariate Normal (MVN) prior is assumed
#' for the columns of \eqn{W}. The multivariate normals have the  0-vector as mean, and diagonal
#' precision matrix \code{Tau}. Different specifications of \code{Tau} (either
#' fixed, or random with Gamma priors) allow achieving regularization on the elements
#' of \eqn{W}. Furthermore, \code{Tau} can be updated with local information, or by sharing
#' information with other elements of the same components of the matrix \eqn{W} (\code{global.var = TRUE}).
#' The latter option can be useful when deciding how many PCs should be used.
#' A fixed \code{Tau} can be updated via Type-II Maximum Likelihood (\code{updatetau=TRUE}).
#' Gamma priors can be activated by setting their scale (\code{alphatau}) and their shape
#' \code{betatau} hyperparameters with values larger than 0. Both \code{alphatau} and \code{betatau}
#' can be scalars (in which case their values are shared across oll PCs), or arrays with component-specific
#' values of these hyperparameters.
#'
#'
#' @param X array_like; \cr
#'          a real \eqn{(I, J)} data matrix (or data frame) to be reduced.
#'
#' @param D integer; \cr
#'          the number of components to be computed.
#'
#' @param nstart integer; \cr
#'               number of (sequential) estimations of the \code{vbpca} algorithm,
#'               with differing random starts.
#'
#' @param maxIter integer; \cr
#'                maximum number of variational algorithm iterations.
#'
#' @param tolerance float; \cr
#'                  stopping criterion for the variational algorithm
#'                  (relative differences between ELBO values).
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
#' @param verbose bool; \cr
#'                logical value which indicates whether the estimation process
#'                information should be printed.
#'
#' @param normalise bool; \cr
#'                  logical argument indicating whether the elements of the weight matrix \eqn{W} should be
#'                  normalised (after model estimation).
#'
#' @param seed integer; \cr
#'             seed used for the random initialization of the model (if \code{svdStart=FALSE} or \code{nstart>1}).
#'
#' @param tau float; \cr
#'            when the prior precisions are fixed, this parameter represents the
#'            values to be used for these priors. When the prior of the precisions are Ingerse-Gamma's,
#'            or \code{updatetau = TRUE}, \code{tau} is the starting value of the precisions.
#'
#' @param updatetau bool; \cr
#'                  when the prior precisions are fixed quantities, this argument specifies whether the
#'                  elements in the precision prior matrix Tau should be updated via Type-II maximum likelihood.
#'
#' @param alphatau float or array_like; \cr
#'                 shape parameter for the Gamma hyperpriors. It can be scalar, or \eqn{D} dimensional.
#'                 The Gamma prior is activated whene \code{alphatau>0} (\code{betatau} must also be
#'                 larger than 0.
#'
#' @param betatau float or array_like; \cr
#'                scale parameter for the Gamma hyperprior. It can be scalar, or \eqn{D} dimensional.
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
#' @param global.var bool; \cr
#'                   it specifies whether \code{tau} should be updated globally (component-specific
#'                   updates) or locally (element-specific updates).
#'
#' @param suppressWarnings bool; \cr
#'                         boolean argument which hides function warnings when set to \code{TRUE}.
#'
#' @param x,object vbpca oject; \cr
#'                  an object of class \code{vbpca}, used as arguments for the \code{print}, \code{is.bayespca} and
#'                  \code{summary} functions.
#'
#' @param ... not used.
#'
#' @return a \code{vbpca} returns a `vbpca` object, which is a list containing the following elements:
#' \item{muW}{ array_like; \cr
#'        posterior means of the weight matrix; \eqn{(J, D)} dimensional array.
#' }
#'
#' \item{P}{ array_like; \cr
#'      point estimate of the (orthogonal) loading matrix; \eqn{(J, D)} dimensional array.
#' }
#'
#' \item{Tau}{ array_like; \cr
#'        the estimates of the prior precisions; depending on the
#'        values of \code{global.var}, it can be a D-dimensional vector or a \eqn{(J, D)} dimensional array.
#' }
#'
#' \item{sigma2}{ float; \cr
#'       point estimate of the variance of the residuals.
#' }
#'
#' \item{HPDI}{ list; \cr
#'       a list containing the high posterior density intervals of the elements of \eqn{W}.
#' }
#'
#' \item{priorAlpha}{ array_like; \cr
#'         a D-dimensional array (or a scalar) containing the values used for the shape hyperparameters of the
#'         Gamma priors.
#' }
#'
#' \item{priorBeta}{ array_like; \cr
#'       a \eqn{(J, D)} or \eqn{D} dimensional array (or a scalar), with the values used for the scale hyperparameters
#'       of the Gamma priors.
#' }
#'
#' \item{elbo}{ float; \cr
#'       evidence lower bound of the model.
#' }
#'
#' \item{converged}{ bool; \cr
#'       boolean denoting whether the Variational Bayes algorithm converged within the required number of iterations.
#' }
#'
#' \item{time}{ array_like; \cr
#'       computation time of the algorithm.
#' }
#'
#' \item{priorvar}{ character; \cr
#'       type of prior for the precisions (either fixed or Gamma).
#' }
#'
#' \item{global.var}{ bool; \cr
#'       \code{global.var} specified as input by the user.
#' }
#'
#' \item{plot}{
#'       traceplot of the evidence lower bounds computed across the various iterations of the algorithm.
#' }
#'
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
#' @author D. Vidotto <d.vidotto@uvt.nl>
#'
#'
#'
#' @examples
#'
#' # Create a synthetic dataset
#' I <- 1e+3
#' X1 <- rnorm(I, 0, 50)
#' X2 <- rnorm(I, 0, 30)
#' X3 <- rnorm(I, 0, 10)
#'
#' X <- cbind(X1, X1, X1, X2, X2, X2, X3, X3 )
#' X <- X + matrix(rnorm(length(X), 0, 1), ncol = ncol(X), nrow = I )
#'
#' # Estimate the Bayesian PCA model, with Gamma priors for tau
#' mod <- vbpca(X, D = 3, alphatau=1e-3, betatau=1e-3 )
#' summary(mod)
#' mod
#'
#'
#'
#'
#'

#' @export
vbpca <- function(X, D = 1, nstart = 1, maxIter = 500, tolerance = 1e-05, center = TRUE, scalecorrection = 1, svdStart = TRUE, verbose = FALSE, normalise = FALSE, seed = 1, 
    tau = 1, updatetau = FALSE, alphatau = 0, betatau = 0, plot.lowerbound = TRUE, hpdi = FALSE, probHPDI = 0.9, global.var = FALSE, suppressWarnings = FALSE) {
    
    UseMethod("vbpca")
    
}

#' @export
is.vbpca <- function(object) {
    class(object) == "vbpca"
}

#' @export
print.vbpca <- function(x, ...) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Final ELBO: ", x$elbo, "\n")
    cat("Elpsed time: ", x$time[[3]], "\n")
}

#' @export
summary.vbpca <- function(object, ...) {
    if (is.null(object) || class(object) != "vbpca") {
        stop("<object> must be a vbpca object.")
    }
    retSum <- c(list(call = object$call), object[1:12])
    class(retSum) <- "summary.vbpca"
    retSum
}


#' @export
print.summary.vbpca <- function(x, ...) {
    cat("Converged: ", x$converged, "\n")
    cat("ELBO: ", x$elbo, "\n")
    cat("Type Prior: ", x$priorvar, "\n")
    cat("Global Variance: ", x$global.var, "\n")
    cat("\n")
    cat("Elapsed time: ", x$time[[3]], "\n")
    
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
}
