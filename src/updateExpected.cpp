#include <RcppArmadillo.h>
#include "aux_functions.h"
#include "expectedvalues.h"

// [[Rcpp::depends(RcppArmadillo)]]
void updateExpected( arma::mat &muW, arma::mat &W2, double &sigma2, arma::mat &muP, int D, 
					 int J, int I, double &hW, const arma::mat& X, const arma::mat& XTX, arma::mat &invTau, bool hpdi, 
					 Rcpp::List &hpdis, double qz, bool scaleprior, double &EWtauW, double v0, 
					 arma::mat &incProbs, double SVS, double &denomX, 
					 double aPostSigma, double normX, 
					 std::string priorvar,
					 bool updatetau, bool globalvar, arma::vec gammatau, arma::vec alphatau,
					 arma::vec &betatau, std::string hypertype, arma::mat &deltataupost,
					 arma::vec deltatau, int JD, arma::vec &priorInclusion, 
					 arma::vec &betastar1, arma::vec &betastar2, bool commonpi, 
					 arma::vec beta1pi, arma::vec beta2pi 
					){

	updateWPXsigma( muW, W2, sigma2, muP, D, J, I, hW, X, XTX,invTau, 
                     priorvar,hpdi,hpdis,qz, scaleprior, EWtauW, v0, 
                     incProbs, SVS, denomX, aPostSigma, normX ); 

	if( priorvar != "fixed" || ( priorvar == "fixed" && updatetau ) ){
		updateTau( sigma2, scaleprior, priorvar, updatetau, SVS, D, invTau, W2, globalvar,
					incProbs, v0, gammatau, alphatau, betatau, hypertype, deltataupost,
					deltatau, JD, J );
	}
	
	if( SVS ){
		updateSVS( incProbs, priorInclusion, betastar1, betastar2, J, D, W2, v0, invTau, sigma2, 
					scaleprior, commonpi, beta1pi, beta2pi, JD );
	}
}
