#ifndef ELBO_FUNCTIONS_H
#define ELBO_FUNCTIONS_H


// Data log-likelihood
void loglikelihood( double &loglik, double sigma2, int J, int I, double denomX );

// Prior for W
void logPriorW( double &logW, arma::mat Tau, int JD, arma::mat W2, int J,
                int D, arma::mat logvar, double invsigma, bool SVS,
                double v0, arma::mat incProbs );

// Prior/Entropy tau
void priorentropyTau( double &logtau, double &hTau, double &priorb, double &hb,
					  int J, int D, arma::mat logvar, bool globalvar,
                      std::string priorvar, arma::mat f, arma::mat Tau,
                      arma::vec alphatau, arma::vec betatau, arma::vec gammatau,
                      arma::vec deltatau, arma::mat deltataupost, std::string hypertype, int JD );
// SVS
void priorentropySVS( double &logPriorIncProbs, double &hPriorIncProbs, double &logPriorGlobalProb,
                      double &hPriorGlobalProb, bool commonpi, arma::vec beta1pi,
                      arma::vec beta2pi, arma::vec priorInclusion, arma::vec betastar1,
                      arma::vec betastar2, arma::mat incProbs, int D );

#endif
