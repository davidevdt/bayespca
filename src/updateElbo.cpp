#include <RcppArmadillo.h>
#include "aux_functions.h"
#include "elbo_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]
void updateElbo( double &elbo, double sigma2, bool globalvar, std::string priorvar,
                 arma::mat W2, arma::vec betatau, arma::vec alphatau,
                 int J, int D,
                 arma::mat Tau, int I, double denomX, double hW, int JD ){

	elbo = 0.0;
	double loglik = 0.0;
	double logW = 0.0;
	double logtau = 0.0;
	double hTau = 0.0;

  arma::mat f = fMat( globalvar, priorvar, W2, 
                      betatau, J, D, alphatau, JD );
  arma::mat logvar = logvarMat( globalvar, J, D, f, priorvar,
   					 			Tau, alphatau );

  loglikelihood( loglik, sigma2, J, I, denomX );

  logPriorW( logW, Tau, JD, W2, J, D, logvar);

	// Note: Entropy for W calculated during E-step

	if( priorvar != "fixed" ){
    priorentropyTau( logtau, hTau, J, D, 
					 logvar, globalvar,
                     f, Tau,
                     alphatau, betatau, JD );
	}
	
	elbo = loglik + logW + logtau + hW + hTau;

}
