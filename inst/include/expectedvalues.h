#ifndef EXPECTEDVALUES_H
#define EXPECTEDVALUES_H

// [[Rcpp::depends(RcppArmadillo)]]

// Xstar, W, P, sigma2
void updateWPXsigma( arma::mat &muW, arma::mat &W2, double &sigma2, arma::mat &muP, int D,
                     int J, int I, double &hW, const arma::mat& X, const arma::mat& XTX, arma::mat Tau,
                     std::string priorvar, bool hpdi, Rcpp::List &hpdis, double qz,
                     double &denomX, double aPostSigma, double normX );

// tau
void updateTau( std::string priorvar, bool updatetau,
                int D, arma::mat &Tau, arma::mat W2, bool globalvar,
                arma::vec alphatau, arma::vec betatau, int JD, int J );

#endif
