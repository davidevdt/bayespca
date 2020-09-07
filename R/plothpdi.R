#' @title Plot high posterior density intervals of vbpca objects. 
#'
#' @aliases plothpdi 
#' 
#' @description Graphical visualization of the HPDI's of the weight matrix computed during the estimation of 
#'              \code{vbpca}. Note that no intervals will be plotted if \code{hpdi = FALSE} is specified in the 
#'              \code{control} argument.  
#'
#' @param obj list; \cr 
#'            a \code{vbpca} object. 
#'
#' @param d integer; \cr 
#'          component number for which the intervals must be plotted. 
#'
#' @param vars array_like; \cr 
#'             an array containing the variable numbers (column numbers) for which the intervals must be plotted. 
#'             If nothing is specified, the method will attempt to plot intervals for all the variables of the data matrix. 
#'
#' @return the graphical visualization of the high posterior density intervals. 
#'
#' @author D. Vidotto <d.vidotto@@uvt.nl>
#' 
#' @seealso \code{\link{vbpca}}
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
#' colnames(X) <- paste('X', 1:8, sep='')
#' # Estimate the Bayesian PCA model, with Inverse Gamma priors for tau
#' # and SVS with Beta priors for priorInclusion, and compute 90% HPD intervals
#' ctrl <- vbpca_control( alphatau = 1., betatau = 1e-02, beta1pi = 1., beta2pi = 1.,
#'                        hpdi = TRUE, probHPDI = 0.9, plot.lowerbound = FALSE  )
#' mod <- vbpca(X, D = 3, priorvar = 'invgamma', SVS = TRUE, control = ctrl )
#' 
#' # Plot the intervals of variables (1, 2, 3) of the second column of the W matrix: 
#' plothpdi(mod, d = 2, vars = 1:3)
#'
#'

#' @export 
plothpdi <- function(obj, d = 1, vars = NULL) {
    
    if (class(obj) != "vbpca") {
        stop("<obj> must be a vbpca object.")
    }
    
    w <- NULL
    
    if (is.null(vars)) {
        vars <- 1:nrow(obj[[1]])
    }
    
    
    if (is.null(obj[[5]])) {
        stop("<hpdi> set to FALSE.")
    }
    
    
    nms <- rownames(obj[[1]])[vars]
    
    dd = data.frame(Variable = nms, w = apply(obj[[5]][[d]][vars, ], 1, mean), lower = obj[[5]][[d]][vars, 1], upper = obj[[5]][[d]][vars, 2])
    
    p <- ggplot(dd, aes(x = dd[, 1], y = w, ymin = dd[, 3], ymax = dd[, 4])) + geom_pointrange() + geom_hline(yintercept = 0, linetype = 2) + coord_flip() + xlab("Variable")
    return(p)
    
    
}
