#ifndef ELBO_FUNCTIONS_H
#define ELBO_FUNCTIONS_H


// Data log-likelihood
void loglikelihood( double &loglik, double sigma2, int J, int I, double denomX );

// Prior for W
void logPriorW( double &logW, arma::mat Tau, int JD, arma::mat W2, int J,
                int D, arma::mat logvar );

// Prior/Entropy tau
void priorentropyTau( double &logtau, double &hTau, 
					  int J, int D, arma::mat logvar, bool globalvar,
                      arma::mat f, arma::mat Tau,
                      arma::vec alphatau, arma::vec betatau, int JD );

#endif
