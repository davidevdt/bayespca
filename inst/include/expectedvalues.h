#ifndef EXPECTEDVALUES_H
#define EXPECTEDVALUES_H

// [[Rcpp::depends(RcppArmadillo)]]

// Xstar, W, P, sigma2
void updateWPXsigma( arma::mat &muW, arma::mat &W2, double &sigma2, arma::mat &muP, int D, 
                     int J, int I, double &hW, const arma::mat& X, const arma::mat& XTX, arma::mat invTau, 
                     std::string priorvar, bool hpdi, Rcpp::List &hpdis, double qz, 
                     bool scaleprior, double &EWtauW, double v0, 
                     arma::mat incProbs, double SVS, double &denomX, 
                     double aPostSigma, double normX );

// tau 
void updateTau( double sigma2, bool scaleprior, std::string priorvar, bool updatetau,
                bool SVS, int D, arma::mat &invTau, arma::mat W2, bool globalvar,
                arma::mat incProbs, double v0, arma::vec gammatau, arma::vec alphatau,
                arma::vec &betatau, std::string hypertype, arma::mat &deltataupost,
                arma::vec deltatau, int JD, int J
               );


// SVS
void updateSVS( arma::mat &incProbs, arma::vec &priorInclusion, arma::vec &betastar1, 
                arma::vec &betastar2, int J, int D, arma::mat W2, double v0, 
                arma::mat invTau, double sigma2, bool scaleprior, 
                bool commonpi, arma::vec beta1pi, arma::vec beta2pi, int JD ); 

#endif
