#' @export
vbpca.default <- function(X, D = 1, nstart = 1, maxIter = 500, tolerance = 1e-05, center = TRUE, scalecorrection = 1, svdStart = TRUE, verbose = FALSE, normalise = FALSE, 
    seed = 1, tau = 1, updatetau = FALSE, alphatau = 0, betatau = 0, plot.lowerbound = TRUE, hpdi = FALSE, probHPDI = 0.9, global.var = FALSE, suppressWarnings = FALSE) {
    
    
    
    ##################################### Convert all parameters
    X <- as.matrix(X)
    D <- as.integer(D)
    nstart <- as.integer(nstart)
    maxIter <- as.integer(maxIter)
    tolerance <- as.double(tolerance)
    center <- as.logical(center)
    scalecorrection <- as.integer(scalecorrection)
    svdStart <- as.logical(svdStart)
    verbose <- as.logical(verbose)
    normalise <- as.logical(normalise)
    seed <- as.double(seed)
    tau <- as.double(tau)
    updatetau <- as.logical(updatetau)
    # priorvar <- as.character(priorvar)
    alphatau <- as.double(alphatau)
    betatau <- as.double(betatau)
    plot.lowerbound <- as.logical(plot.lowerbound)
    hpdi <- as.logical(hpdi)
    probHPDI <- as.double(probHPDI)
    global.var <- as.logical(global.var)
    
    
    ##################################### Checks
    if (D <= 0 || D > ncol(X)) {
        stop("The number of components <D> must be > 0 and <= ncol(X).")
    }
    
    if (tolerance <= 0) {
        stop("The convergence criterion <tolerance> must be > 0.")
    }
    
    if (maxIter <= 0) {
        stop("The max number of iterations <maxIter> must be > 0.")
    }
    
    if (nstart <= 0) {
        stop("The number of random starts <nstart> must be > 0.")
    }
    
    if (length(tau) > 1) {
        stop("<tau>  must be scalar.")
    }
    
    if (tau <= 0) {
        stop("<tau> must be > 0.")
    }
    
    if (all(alphatau > 0) && all(betatau > 0)) {
        priorvar <- "gamma"
        updatetau <- TRUE 
    } else {
        priorvar <- "fixed"
        alphatau <- 0
        betatau <- 0
    }
    
    
    if (global.var == FALSE) {
        
        if (priorvar == "gamma") {
            
            if (length(alphatau) == 1) {
                alphatau <- rep(alphatau, D)
            }
            
            if (length(betatau) == 1) {
                betatau <- rep(betatau, D)
            }
            
            if (any(c(length(alphatau), length(betatau)) != D)) {
                stop("The size of <alphatau> and <betatau> must be either 1 or D.")
            }
            
        }
        
        # Method information
        if (verbose) {
            
            if (priorvar == "fixed") {
                
                if (!updatetau) {
                  message("Local prior variances : fixed.")
                } else {
                  message("Local prior variances : fixed, Type-II ML update.")
                }
                
            } else {
                
                message("Local prior variances : Gamma.")
            }
        }
        
    } else {
        
        if (priorvar == "gamma") {
            
            if (length(alphatau) == 1) {
                alphatau <- rep(alphatau, D)
            }
            
            if (length(betatau) == 1) {
                betatau <- rep(betatau, D)
            }
            
            if (any(c(length(alphatau), length(betatau)) != D)) {
                stop("The size of <alphatau> and <betatau> must be either 1 or D.")
            }
            
        }
        
        # Method information
        if (verbose) {
            if (priorvar == "fixed") {
                
                if (!updatetau) {
                  message("Local prior variances : fixed.")
                } else {
                  message("Local prior variances : fixed, Type-II ML update.")
                }
                
            } else {
                
                message("Global prior variances: Gamma.")
                
            }
        }
    }
    
    
    ##################################### HPD intervals
    if (probHPDI <= 0 | probHPDI >= 1) {
        stop("<probHPDI> must be in the interval (0,1).")
    }
    
    # Quantiles
    qz <- qnorm(1 - ((1 - probHPDI)/2))
    
    
    
    ##################################### Initializations and scaling
    X <- na.omit(X)
    J <- ncol(X)
    I <- nrow(X)
    
    if (center) {
        means <- apply(X, 2, mean)
        X <- t(t(X) - means)
        rm(means)
    }
    
    if (scalecorrection >= 0) {
        
        sdX <- apply(X, 2, function(y) sqrt(sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)/(I - scalecorrection)))
        X <- t(t(X)/sdX)
        rm(sdX)
        
    }
    
    # Precision Matrix
    JD <- J * D
    Tau <- matrix(tau, J, D)
    
    
    if (seed > 0) {
        set.seed(seed)
    }
    
    
    ##################################### Model Estimation
    timeStart <- proc.time()
    
    retList <- mainBayesPCA(X, D, I, J, nstart, maxIter, tolerance, svdStart, verbose, updatetau, priorvar, alphatau, betatau, JD, Tau, qz, global.var, hpdi)
    
    finalTime <- proc.time() - timeStart
    rm(timeStart)
    retList$time <- finalTime
    
    
    ##################################### Results evaluation
    if (!suppressWarnings) {
        try(if (!retList$globalConverged) 
            warning("vbpca has not converged. Please re-run by increasing <maxIter> or the convergence criterion <tolerance>.", call. = FALSE, immediate. = FALSE, noBreaks. = FALSE, 
                domain = NULL))
        
        if (retList$globalConverged) {
            
            MX <- max(retList$elbovals[-1])
            try(if (retList$elbovals[length(retList$elbovals)] != MX) 
                warning("decreasing ELBO.", call. = FALSE, immediate. = FALSE, noBreaks. = FALSE, domain = NULL))
            rm(MX)
        }
    }
    
    
    if (normalise) {
        retList$globalMuW = apply(retList$globalMuW, 2, function(x) x/sqrt(sum(x^2)))
    }
    
    
    if (!is.null(colnames(X))) {
        nms <- colnames(X)
    } else {
        nms <- paste("variable", 1:J)
    }
    
    
    if (!hpdi) {
        
        retList$globalHPDIS <- NULL
        
    } else {
        
        if (retList$globalConverged) {
            names(retList$globalHPDIS) <- noquote(paste("Component ", 1:D, sep = ""))
            
            for (i in seq_along(retList$globalHPDIS)) {
                colnames(retList$globalHPDIS[[i]]) <- noquote(paste(c(((100 - (probHPDI * 100))/2), 100 - ((100 - probHPDI * 100)/2)), "%", sep = ""))
                rownames(retList$globalHPDIS[[i]]) <- nms
            }
            
        }
        
    }
    
    
    if (global.var) {
        
        retList$globaltau <- retList$globaltau[1, ]
        
    } else {
        
        colnames(retList$globaltau) <- paste("Component", 1:D)
        rownames(retList$globaltau) <- nms
        
    }
    rm(tau)
    
    
    pl <- NULL
    
    if (retList$globalConverged) {
        if (plot.lowerbound & (global.var)) {
            
            par(mfrow = c(2, 1))
            plot(retList$elbovals[-1], type = "l", col = "blue", lwd = 2, main = "Evidence Lower Bound", ylab = "ELBO", xlab = "Iteration")
            plot(retList$globaltau, type = "b", col = "blue", lwd = 2, main = "Prior Precisions", ylab = expression(tau), xlab = "Component")
            
            pl <- recordPlot()
            par(mfrow = c(1, 1))
            
            
        } else if (plot.lowerbound & !(global.var)) {
            
            plot(retList$elbovals[-1], type = "l", col = "blue", lwd = 2, main = "Evidence Lower Bound", ylab = "ELBO", xlab = "Iteration")
            
            pl <- recordPlot()
            
        } else if (!plot.lowerbound & (global.var)) {
            
            plot(retList$globaltau, type = "b", col = "blue", lwd = 2, main = "Prior Precisions", ylab = expression(tau), xlab = "Component")
            
            pl <- recordPlot()
            par(mfrow = c(1, 1))
        }
    }
    
    colnames(retList$globalMuW) <- paste("Component", 1:D)
    colnames(retList$globalMuP) <- paste("Component", 1:D)
    rownames(retList$globalMuW) <- nms
    rownames(retList$globalMuP) <- nms
    
    
    ##################################### Output
    ret <- list(muW = retList$globalMuW, P = retList$globalMuP, Tau = retList$globaltau, sigma2 = retList$globalSigma2, HPDI = retList$globalHPDIS, priorAlpha = alphatau, 
        priorBeta = betatau, elbo = retList$globalElbo, converged = retList$globalConverged, time = retList$time, priorvar = priorvar, global.var = global.var, plot = invisible(pl))
    ret$call <- match.call()
    class(ret) <- "vbpca"
    ret
}
