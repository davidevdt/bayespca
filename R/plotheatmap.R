#' @title Plot weights/prior precision matrix heatmaps.
#'
#' @aliases plotheatmap
#'
#' @description Heatmap of the estimated weights (W) or prior precision matrix (Tau) estimated with \code{vbpca}.
#'              Furthermore, the function allows specifying a truncating (usually large) value for the prior precision
#'              matrix, so that all elements in the matrix exeeding such value will be truncated. This will also cause the
#'              corresponding elements of the weights matrix to be set to 0, making the matrix sparse. The function will
#'              also return such matrices.
#'
#' @param obj list; \cr
#'            a \code{vbpca} object.
#'
#' @param matrix_type character; \cr
#'                    a character vector defining the matrix that the user desires to plot; either 'W'
#'                    (for the weights matrix) or 'Tau' (for the prior precision matrix)
#'
#' @param bound_tau float; \cr
#'                  truncating value for the prior precision matrix. This will not only cause truncation in the
#'                  Tau matrix, but also will make the weights matrix W sparse.
#'
#' @param limits_heatmap array_like; \cr
#'                       2-dimensional array containing limits for the scale_fill_gradient2 functionality of ggplo2.
#'                       Meaningful only when \code{matrix_type = 'W'}.
#'
#' @return the heatmap of the selected matrix, as well as the (truncated) prior precision matrix Tau, and the (sparse)
#'          weights matrix W.
#'
#' @author D. Vidotto <d.vidotto@uvt.nl>
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
#' # Estimate the Bayesian PCA model, with Gamma priors for tau
#' mod <- vbpca(X, D = 3, alphatau=1e-3, betatau=1e-3 )
#'
#' # Plot heatmap of Tau (setting a threshold value at 50):
#' plotheatmap(mod, matrix_type='Tau', bound_tau=50)
#'
#'

#' @export
plotheatmap <- function(obj, matrix_type = "Tau", bound_tau = NULL, limits_heatmap = NULL) {
    
    if (class(obj) != "vbpca") {
        stop("<obj> must be a vbpca object.")
    }
    
    if (matrix_type != "Tau" && matrix_type != "W") {
        stop("<matrix_type> must be either 'Tau' (for heatmap of precision matrix) or 'W' (for heatmap of weights).")
    }
    
    col_nms <- colnames(obj[[1]])
    row_nms <- rownames(obj[[1]])
    
    ret_W <- obj[[1]]
    ret_Tau <- obj[[3]]
    
    if (!is.null(bound_tau)) {
        
        ret_Tau[obj[[3]] > bound_tau] <- bound_tau
        ret_W[obj[[3]] > bound_tau] <- 0
        
    }
    
    if (matrix_type == "W") {
        y_ax <- x_ax <- vals <- NULL 
        dd <- expand.grid(y_ax = row_nms, x_ax = col_nms)
        dd$vals <- c(ret_W)
        
        if (is.null(limits_heatmap)) {
            
            limits_heatmap <- c(min(dd$vals), max(dd$vals))
            
        } else {
            
            if (length(limits_heatmap) != 2) {
                stop("<limits_heatmap> must be of length 2.")
            }
            
        }
        
        p <- ggplot(dd, aes(x_ax, ordered(y_ax, levels = rev(sort(unique(y_ax)))), fill = vals)) + geom_tile() + scale_fill_gradient2(low = "red", mid = "black", high = "green", 
            midpoint = 0, limits = limits_heatmap, oob = squish) + ylab("") + xlab("") + ggtitle("W Matrix") + theme_bw() + theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5), 
            panel.grid.major = element_blank(), panel.border = element_blank())
    } else {
        
        dd <- expand.grid(y_ax = row_nms, x_ax = col_nms)
        dd$vals <- c(ret_Tau)
        
        p <- ggplot(dd, aes(x_ax, ordered(y_ax, levels = rev(sort(unique(y_ax)))), fill = vals)) + geom_tile() + scale_fill_distiller(palette = "Blues", direction = +1) + 
            ylab("") + xlab("") + ggtitle("Tau Matrix") + theme_bw() + theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), 
            panel.border = element_blank())
        
    }
    
	print(p) 
	
    return(invisible(list(W_matrix = ret_W, Tau_matrix = ret_Tau)))
    
}
