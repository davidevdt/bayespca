#include <RcppArmadillo.h>
#include "updateExpected.h"
#include "updateElbo.h"
#include "aux_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]
void vbalgorithm( const arma::mat& X, int D, int I, int J, int maxIter, double tolerance,
                  bool updatetau, std::string priorvar, arma::vec alphatau, arma::vec betatau,
                  int JD, arma::mat Tau, double qz, Rcpp::List hpdis, int &it, double sigma2,
                  double aPostSigma, double denomX, double normX, double elbo,
				          double hW, double &finalElbo, bool &converged,
                  arma::mat muW, arma::mat muP, arma::mat W2, const arma::mat& XTX,
                  double &globalElbo, arma::mat &globalMuW, arma::mat &globalMuP,
                  arma::mat &globaltau, double &globalSigma2,
                  bool &globalconverged, arma::vec &elbovals,
                  Rcpp::List &globalhpdi, bool globalvar, bool hpdi, bool verbose ){

	finalElbo = 0.0;
	arma::vec allElbos;
  arma::vec elboValue = arma::vec( maxIter + 1);

	elboValue.fill( -(arma::datum::inf) );
	converged = FALSE;

	for( it = 1; it < (maxIter + 1); it++ ){

    updateExpected( muW, W2, sigma2, muP, D, J, I, hW, X, XTX,
    				Tau, hpdi, hpdis, qz, denomX,
    				aPostSigma, normX, priorvar,
    				updatetau, globalvar, JD,  
					alphatau, betatau );


    updateElbo( elbo, sigma2, globalvar, priorvar,
                W2, betatau, alphatau, J, D,
                Tau, I, denomX, hW, JD );

		elboValue( it ) = elbo;

		if( std::abs( (elboValue( it ) - elboValue( it - 1 ))/elboValue( it ) ) <= tolerance ){
			converged = TRUE;
			break;
		}

		if( it % 10 == 0 ){
			Rcpp::checkUserInterrupt();
		}

		if( verbose & ((it - 1) % 10 == 0) ){
			Rcpp::Rcout << "Iteration: " << it << " - ELBO: " << elboValue(it) << std::endl;
		}
	}

	allElbos = elboValue( find( elboValue != -(arma::datum::inf) ) );
	finalElbo = allElbos( allElbos.size() - 1 );

	if( finalElbo > globalElbo ){
		globalElbo = finalElbo;
		globalMuW = muW;
		globalMuP = muP;
		globaltau = Tau;
		globalSigma2 = sigma2;
		globalconverged = converged;
		globalhpdi = hpdis;
		elbovals = allElbos;
	}
}
