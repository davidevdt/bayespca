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
#' @param X array_like; \cr 
#'          the data matrix used to estimate the \code{vbpca} object. 
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
#' \dontrun{
#'
#' # mod is a vbpca object estimated on the X matrix; 
#' # to plot the intervals of elements (1, 2, 3) of the second column of the W matrix: 
#' plothpdi(mod, d = 2, X, vars = 1:3)
#' }
#'
#'

#' @export 
plothpdi <- function(obj, d = 1, X, vars = NULL) {
    
    if (class(obj) != "vbpca") {
        stop("<obj> must be a vbpca object.")
    }
    
    w <- NULL
    
    if (is.null(vars)) {
        vars <- 1:ncol(X)
    }
    
    
    if (is.null(obj[[5]])) {
        stop("<hpdi> set to FALSE.")
    }
    
    X <- X[, vars]
    
    if (is.null(colnames(X))) {
        colnames(X) = paste("variable", vars)
    }
    
    
    dd = data.frame(Variable = colnames(X), w = apply(obj[[5]][[d]][vars, ], 1, mean), lower = obj[[5]][[d]][vars, 1], upper = obj[[5]][[d]][vars, 2])
    
    p <- ggplot(dd, aes(x = dd[, 1], y = w, ymin = dd[, 3], ymax = dd[, 4])) + geom_pointrange() + geom_hline(yintercept = 0, linetype = 2) + coord_flip() + xlab("Variable")
    return(p)
    
    
}
