#' @export 
vbpca.default <- function(X, D = 1, maxIter = 500, tolerance = 1e-05, verbose = FALSE, tau = 1, updatetau = FALSE, priorvar = "invgamma", SVS = FALSE, priorInclusion = 0.5, global.var = FALSE, control = list(), 
    suppressWarnings = FALSE) {
    
    
    
    ##################################### Load and read all parameters
    
    ### Control parameters
    ctrl <- vbpca_control()
    namesCtrl <- names(ctrl)
    
    ctrl[(namesSpec <- names(control))] <- control
    
    if (length(unkNames <- namesSpec[!namesSpec %in% namesCtrl])) {
        warning("unknown names in control: ", paste(unkNames, collapse = ", "))
    }
    
    nstart <- ctrl$nstart
    center <- ctrl$center
    scalecorrection <- ctrl$scalecorrection
    svdStart <- ctrl$svdStart
    normalise <- ctrl$normalise
    seed <- ctrl$seed
    plot.lowerbound <- ctrl$plot.lowerbound
    hpdi <- ctrl$hpdi
    probHPDI <- ctrl$probHPDI
    scaleprior <- ctrl$scaleprior
    alphatau <- ctrl$alphatau
    betatau <- ctrl$betatau
    gammatau <- ctrl$gammatau
    deltatau <- ctrl$deltatau
    hypertype <- ctrl$hypertype
    beta1pi <- ctrl$beta1pi
    beta2pi <- ctrl$beta2pi
    v0 <- ctrl$v0
    ### 
    
    
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
    priorvar <- as.character(priorvar)
    alphatau <- as.double(alphatau)
    betatau <- as.double(betatau)
    gammatau <- as.double(gammatau)
    deltatau <- as.double(deltatau)
    SVS <- as.logical(SVS)
    priorInclusion <- as.double(priorInclusion)
    beta1pi <- as.double(beta1pi)
    beta2pi <- as.double(beta2pi)
    v0 <- as.double(v0)
    plot.lowerbound <- as.logical(plot.lowerbound)
    hpdi <- as.logical(hpdi)
    probHPDI <- as.double(probHPDI)
    scaleprior <- as.logical(scaleprior)
    hypertype <- as.character(hypertype)
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
    
    
    if (global.var == FALSE) {
        
        if (priorvar != "fixed" & priorvar != "invgamma" & priorvar != "jeffrey") {
            stop("<priorvar> for local variance priors: either 'fixed', 'invgamma', or 'jeffrey'.")
        }
        
        if (priorvar == "invgamma" & all(gammatau > 0)) {
            if (hypertype != "common" & hypertype != "component" & hypertype != "local") {
                stop("<hypertype>: either 'common', 'component' or 'local'.")
            }
        }
        
        if (priorvar == "invgamma") {
            
            if (any(alphatau <= 0) | any(betatau <= 0)) {
                stop("<alphatau> and <betatau> must be larger than 0.")
            }
            
            
            if (length(alphatau) == 1) {
                alphatau <- rep(alphatau, D)
            }
            
            if (length(betatau) == 1) {
                betatau <- rep(betatau, D)
            }
            
            if (length(alphatau) != 1 && hypertype == "common") {
                alphatau <- rep(alphatau[1], D)
            }
            
            
            if (length(betatau) != 1 && hypertype == "common") {
                betatau <- rep(betatau[1], D)
            }
            
            
            if (any(c(length(alphatau), length(betatau)) != D)) {
                stop("The size of <alphatau> and <betatau> must be either 1 or D.")
            }
            
            
            if (all(gammatau > 0)) {
                
                if (any(deltatau <= 0)) {
                  stop("<deltatau> must be larger than 0.")
                }
                
                if (length(gammatau) == 1) {
                  gammatau <- rep(gammatau, D)
                }
                
                if (length(deltatau) == 1) {
                  deltatau <- rep(deltatau, D)
                }
                
                if (length(gammatau) != 1 & hypertype == "common") {
                  gammatau <- rep(gammatau[1], D)
                }
                
                if (length(deltatau) != 1 & hypertype == "common") {
                  deltatau <- rep(deltatau[1], D)
                }
                
                if (any(c(length(gammatau), length(deltatau)) != D)) {
                  stop("The size of <alphatau>, <betatau>, <gammatau>, <deltatau> must be either 1 or D.")
                }
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
                
            } else if (priorvar == "jeffrey") {
                
                message("Local prior variances : Jeffreys' prior.")
                
            } else {
                
                if (any(gammatau <= 0)) {
                  message("Local prior variances : Inverse-Gamma, fixed hyperparameters.")
                } else {
                  if (hypertype == "common") {
                    message("Local prior variances : Inverse-Gamma, random common hyperparameters.")
                  } else if (hypertype == "component") {
                    message("Local prior variances : Inverse-Gamma, random component-specific hyperparameters.")
                  } else {
                    message("Local prior variances : Inverse-Gamma, random local hyperparameters.")
                  }
                }
                
            }
        }
    } else {
        # global prior.var
        
        if (priorvar != "fixed" & priorvar != "invgamma" & priorvar != "jeffrey") {
            stop("<priorvar> for global variance priors: either 'fixed', 'invgamma', or 'jeffrey'.")
        }
        
        
        if (priorvar == "invgamma") {
            
            hypertype <- "component"
            
            if (any(c(alphatau, betatau) <= 0)) {
                stop("<alphatau> and <betatau> must be >0.")
            }
            
            if (length(alphatau) == 1) {
                alphatau <- rep(alphatau, D)
            }
            if (length(betatau) == 1) {
                betatau <- rep(betatau, D)
            }
            
            if (any(c(length(alphatau), length(betatau)) != D)) {
                stop("The size of <alphatau> and <betatau> must be either 1 or D.")
            }
            
            if (all(gammatau > 0)) {
                
                if (any(deltatau <= 0)) {
                  stop("<deltatau> must be >0.")
                }
                
                if (length(gammatau == 1)) {
                  gammatau <- rep(gammatau, D)
                }
                if (length(deltatau == 1)) {
                  deltatau <- rep(deltatau, D)
                }
                if (any(c(length(gammatau), length(deltatau)) != D)) {
                  stop("The size of <alphatau>, <betatau>, <gammatau>, and <deltatau> must be either 1 or D.")
                }
                
            }
        }
        
        # Method information
        if (verbose) {
            if (priorvar == "fixed") {
                
                message("Global prior variances: fixed.")
                
            } else if (priorvar == "jeffrey") {
                
                message("Global prior variances: Jeffrey's prior.")
                
            } else {
                
                if (any(gammatau <= 0)) {
                  message("Global prior variances: Inverse-Gamma, fixed hyperparameters.")
                } else {
                  message("Global prior variances: Inverse-Gamma, (component-specific) random hyperparameters.")
                }
            }
        }
    }
    
    
    ##################################### Stochastic Variable Selection
    commonpi <- TRUE
    
    if (SVS) {
        
        if ((v0 <= 0) | (v0 >= 1)) {
            stop("<v0> should be in the interval (0,1), preferably close to 0.")
        }
        
        
        if (any(priorInclusion <= 0)) {
            stop("<priorInclusion> must be in the interval (0,1).")
        }
        
        
        
        if (length(priorInclusion) == D) {
            commonpi <- FALSE
        }
        
        if (length(priorInclusion) == 1) {
            priorInclusion <- rep(priorInclusion, D)
        }
        
        # Beta priors on inclusion probabilities
        if (all(beta1pi > 0)) {
            if (any(beta2pi <= 0)) {
                stop("<beta2pi> must be larger than 0.")
            }
            
            if (all(c(length(beta1pi), length(beta2pi)) == D)) {
                commonpi <- FALSE
            }
            if (length(beta1pi) == 1) {
                beta1pi <- rep(beta1pi, D)
            }
            if (length(beta2pi) == 1) {
                beta2pi <- rep(beta2pi, D)
            }
            
            if (length(beta1pi) != 1 & commonpi) {
                beta1pi <- rep(beta1pi[1], D)
            }
            if (length(beta2pi) != 1 & commonpi) {
                beta2pi <- rep(beta2pi[1], D)
            }
            
            if (any(c(length(beta1pi), length(beta2pi)) != D)) {
                stop("The size of <beta1pi> and <beta2pi> must be either 1 or D.")
            }
        }
        
        # Method information
        if (verbose) {
            if (commonpi) {
                if (any(beta1pi <= 0)) {
                  
                  message("SVS activated: common prior probabilities, fixed.")
                  
                } else {
                  
                  message("SVS activated: random common prior probabilities with Beta priors.")
                  
                }
            } else {
                if (any(beta1pi <= 0)) {
                  
                  message("SVS activated: component-specific prior probabilities, fixed.")
                  
                } else {
                  
                  message("SVS activated: random component-specific prior probabilities with Beta priors.")
                  
                }
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
        
        sdX <- apply(X, 2, function(y) sqrt(sum(y^2, na.rm = TRUE)/(I - scalecorrection)))
        X <- t(t(X)/sdX)
        rm(sdX)
        
    } else {
        
        if (!suppressWarnings) {
            warning("unscaled data - ELBO values might be positive.", call. = FALSE, immediate. = TRUE, noBreaks. = FALSE, domain = NULL)
        }
        
    }
    
    # Inverse variances
    JD <- J * D
    invTau <- matrix(1/tau, J, D)
    
    
    if (seed > 0) {
        set.seed(seed)
    }
    
    
    ##################################### Model Estimation
    timeStart <- proc.time()
    
    retList <- mainBayesPCA(X, D, I, J, nstart, maxIter, tolerance, svdStart, verbose, updatetau, priorvar, alphatau, betatau, gammatau, deltatau, SVS, priorInclusion, beta1pi, beta2pi, v0, commonpi, JD, invTau, 
        qz, scaleprior, hypertype, global.var, hpdi)
    
    finalTime <- proc.time() - timeStart
    rm(timeStart)
    retList$time <- finalTime
    
    
    
    ##################################### Results evaluation
    if (!suppressWarnings) {
        try(if (!retList$globalConverged) 
            warning("vbpca has not converged. Please re-run by increasing <maxIter> or the convergence criterion <tolerance>.", call. = FALSE, immediate. = FALSE, noBreaks. = FALSE, domain = NULL))
        
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
            plot(1/retList$globaltau, type = "b", col = "blue", lwd = 2, main = "Prior Variances", ylab = expression(tau), xlab = "Component")
            
            pl <- recordPlot()
            par(mfrow = c(1, 1))
            
            
            
        } else if (plot.lowerbound & !(global.var)) {
            
            plot(retList$elbovals[-1], type = "l", col = "blue", lwd = 2, main = "Evidence Lower Bound", ylab = "ELBO", xlab = "Iteration")
            
            pl <- recordPlot()
            
        } else if (!plot.lowerbound & (global.var)) {
            
            plot(1/retList$globaltau, type = "b", col = "blue", lwd = 2, main = "Prior Variances", ylab = expression(tau), xlab = "Component")
            
            pl <- recordPlot()
        }
        
        
        
    }
    
    
    
    colnames(retList$globalMuW) <- paste("Component", 1:D)
    colnames(retList$globalMuP) <- paste("Component", 1:D)
    rownames(retList$globalMuW) <- nms
    rownames(retList$globalMuP) <- nms
    
    
    
    if (!SVS) {
        retList$globalPriorInc <- NULL
        retList$globalIncPr <- NULL
    } else {
        colnames(retList$globalIncPr) <- paste("Component", 1:D)
        rownames(retList$globalIncPr) <- nms
    }
    
    
    
    if (priorvar == "invgamma") {
        
        if (all(gammatau > 0)) {
            if (hypertype == "common") {
                names(retList$globalbetatau) <- paste("beta", 1:D)
            } else if (hypertype == "component") {
                names(retList$globalbetatau) <- paste("beta", 1:D)
            } else {
                colnames(retList$globalbetatau) <- paste("beta", 1:D)
                rownames(retList$globalbetatau) <- nms
            }
        }
        
    } else {
        
        retList$globalbetatau <- betatau
        
    }
    
    
    ##################################### Output
    ret <- list(muW = retList$globalMuW, P = retList$globalMuP, invTau = retList$globaltau, sigma2 = retList$globalSigma2, HPDI = retList$globalHPDIS, priorAlpha = alphatau, priorBeta = retList$globalbetatau, 
        priorInclusion = retList$globalPriorInc, inclusionProbabilities = retList$globalIncPr, elbo = retList$globalElbo, converged = retList$globalConverged, time = retList$time, priorvar = priorvar, global.var = global.var, 
        hypertype = hypertype, SVS = SVS, plot = invisible(pl))
    
    ret$call <- match.call()
    
    class(ret) <- "vbpca"
    ret
    
}
