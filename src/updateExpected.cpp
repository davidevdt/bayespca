#include <RcppArmadillo.h>
#include "aux_functions.h"
#include "expectedvalues.h"

// [[Rcpp::depends(RcppArmadillo)]]
void updateExpected( arma::mat &muW, arma::mat &W2, double &sigma2, arma::mat &muP, int D,
										 int J, int I, double &hW, const arma::mat& X, const arma::mat& XTX,
										 arma::mat &Tau, bool hpdi, Rcpp::List &hpdis, double qz, double &denomX,
										 double aPostSigma, double normX, std::string priorvar,
										 bool updatetau, bool globalvar, int JD, 
										 arma::vec alphatau, arma::vec betatau ){

	updateWPXsigma( muW, W2, sigma2, muP, D, J, I,
									hW, X, XTX, Tau, priorvar, hpdi, hpdis,
									qz, denomX, aPostSigma, normX );
									
	if( updatetau ){
		updateTau(priorvar, updatetau,
	             D, Tau, W2, globalvar,
	             alphatau, betatau, JD, J );
	}
	
}
