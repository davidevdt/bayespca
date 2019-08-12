mainBayesPCA <- function(X, D, I, J, nstart, maxIter, tolerance, svdStart, verbose, updatetau, priorvar, alphatau, betatau, gammatau, deltatau, SVS, priorInclusion, beta1pi, beta2pi, v0, commonpi, JD, invTau, qz, scaleprior, hypertype, globalvar, hpdi) {
    .Call(`_bayespca_mainBayesPCA`, X, D, I, J, nstart, maxIter, tolerance, svdStart, verbose, updatetau, priorvar, alphatau, betatau, gammatau, deltatau, SVS, priorInclusion, beta1pi, beta2pi, v0, commonpi, JD, invTau, qz, scaleprior, hypertype, globalvar, hpdi)
}

