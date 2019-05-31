#include <RcppArmadillo.h>
#include "aux_functions.h"
#include "elbo_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]
void updateElbo( double &elbo, double sigma2, bool globalvar, std::string priorvar, 
                 arma::mat W2, bool SVS, arma::mat incProbs, double v0, arma::vec gammatau,
                 arma::vec betatau, arma::mat deltataupost, arma::vec alphatau, 
                 arma::vec deltatau, std::string hypertype, int J, int D, 
                 arma::mat invTau, int I, double denomX, double hW, int JD, 
                 bool commonpi, arma::vec beta1pi, arma::vec beta2pi, arma::vec priorInclusion, 
                 arma::vec betastar1, arma::vec betastar2, bool scaleprior ){

	elbo = 0.0; 
	double loglik = 0.0; 
	double logW = 0.0; 
	double logtau = 0.0; 
	double priorb = 0.0; 
	double hb = 0.0; 
	double hTau = 0.0; 
	double logPriorIncProbs = 0.0; 
	double hPriorIncProbs = 0.0; 
	double logPriorGlobalProb = 0.0; 
	double hPriorGlobalProb = 0.0; 
	
	double tmpinvsigma; 
	
	if( scaleprior ){
		tmpinvsigma = 1/sigma2; 
	}else{
	  tmpinvsigma = 1.0; 
	}

	arma::mat f = fMat( globalvar, priorvar, W2, tmpinvsigma, SVS, incProbs, v0, 
						gammatau, betatau, deltataupost, hypertype, 
						J, D, alphatau, JD ); 

	
	arma::mat logvar = logvarMat( globalvar, J, D, f, priorvar, invTau, 
									   hypertype, alphatau );
	
	
	loglikelihood( loglik, sigma2, J, I, denomX );  
	
 
	logPriorW( logW, invTau, JD, W2, J, D, logvar, tmpinvsigma, SVS, 
			   v0, incProbs ); 
			  
			  
	// Entropy for W calculated during E-step
	
	if( priorvar != "fixed" ){
    priorentropyTau( logtau, hTau, priorb, hb, J, D, logvar, globalvar, 
                     priorvar, f, invTau, alphatau, betatau, 
                     gammatau, deltatau, deltataupost, hypertype, JD  );
	}
	

	
	if( SVS ){		
		priorentropySVS( logPriorIncProbs, hPriorIncProbs, logPriorGlobalProb, 
						 hPriorGlobalProb, commonpi, beta1pi, beta2pi, 
						 priorInclusion, betastar1, betastar2,
						 incProbs, D ); 
	}	
	
	elbo = loglik + logW + logtau + priorb + logPriorIncProbs + logPriorGlobalProb + hW + hTau + hb + hPriorIncProbs + hPriorGlobalProb; 
	
}