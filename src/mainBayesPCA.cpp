/*

	mainBayesPCA: estimate Variational Bayes PCA with C++.
	The routine uses the RcppArmadillo library to build and
	compute its linear algebra functions.

	The vbalgorithm() function is called for the implementation
	of the proper algorithm. In turn, vbalgorithm() is mainly
	composed of two steps: (a) updateExpected(), which updates
	the expected values of the latent variables (muW) and
	the point estimates of the other model parameters (muP,
	sigma2, tau, and so on ); (b) updateElbo, which computes
	the evidence lower bound (ELBO) given the current parameter
	estimates.


*/

#include <RcppArmadillo.h>
#include "vbalgorithm.h"
#include "aux_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List mainBayesPCA(	const arma::mat& X, int D, int I, int J, int nstart, int maxIter,
                          double tolerance, bool svdStart, bool verbose,
                          bool updatetau, std::string priorvar,
                          arma::vec alphatau, arma::vec betatau,
                          int JD, arma::mat Tau, double qz, bool globalvar, bool hpdi
              		       ){

	// arma::arma_rng::set_seed( seed );

	// Initializations : Lists
	Rcpp::List ret;
	Rcpp::List hpdis;

	if( hpdi ){
		hpdis = Rcpp::List(D);
	}else{
		hpdis = Rcpp::List(1);
	}

	// Initializations: Integers
	int st;
	int it;

	// Initializations : Doubles
	double sigma2 = 1.0 ;
	double aPostSigma;
	aPostSigma = ( double(J) * double(I));
	double denomX = 0.0;
	double normX;
	double elbo = 0.0;
	double hW = 0.0;
	double finalElbo = 0.0;

	// Initializations: Booleans
	bool converged = FALSE;

	// Initializations : Matrices
	arma::mat muW(J, D);
	arma::mat muP(J, D);
	arma::mat W2(J, D);
	arma::mat XTX(J, J);

	// Calculate XTX and norm(X)
	funcX( X, XTX, normX );

  // Initializations: muW and muP
	if( svdStart ){
		arma::mat U1(I, D);
		arma::vec d1(D);
		SVD( X, U1, d1, muW, D, D );
		//muP = Rcpp::clone(muW) ;
		muP = muW ;
	}else{
		muW = arma::randu(J, D);
		muP = arma::randu(J, D);
	}

  //Initializations: Outputs
	double globalElbo = - ( arma::datum::inf );
	arma::mat globalMuW( J , D, arma::fill::zeros );
	arma::mat globalMuP( J, D, arma::fill::zeros );
	arma::mat globaltau(J, D, arma::fill::zeros);
	double globalSigma2 = 0.0;
	bool globalconverged = FALSE;
	arma::vec elbovals;
	Rcpp::List globalhpdi;

	// Run the Algorithm
	for( st = 0; st < nstart; st++ ){

		// Random starts
		if( st > 0){

			if( svdStart ){
				muW += (arma::randn(J, D) * 0.05);
				muP += (arma::randn(J, D) * 0.05);
			}else{
				muW = arma::randu(J, D);
				muP = arma::randu(J, D);
			}
		}

		W2 = arma::square( muW );

		vbalgorithm(X, D, I, J, maxIter, tolerance, updatetau,
      					priorvar, alphatau, betatau, JD, Tau, qz, hpdis,
      					it, sigma2, aPostSigma, denomX, normX, elbo, hW,
      					finalElbo, converged, muW, muP, W2, XTX,
                globalElbo, globalMuW, globalMuP, globaltau,
      					globalSigma2, globalconverged, elbovals,
      					globalhpdi, globalvar, hpdi, verbose );


		if( verbose ){
			if( converged ){
				Rcpp::Rcout << "Start # " << st + 1 << " has converged in " << (it-1) << " iterations; lower bound = " << finalElbo << std::endl;
			}else{
				Rcpp::Rcout << "Start # " << st + 1 << " has not converged after " << (it-1) << " iterations." << std::endl;
			}
		}

	}

	ret["globalElbo"] = globalElbo;
	ret["globalMuW"] = globalMuW;
	ret["globalMuP"] = globalMuP;
	ret["globaltau"] = globaltau;
	ret["globalSigma2"] = globalSigma2;
	ret["globalConverged"] = globalconverged;
	ret["globalHPDIS"] = globalhpdi;
	ret["elbovals"] = elbovals;

	return ret;
}
