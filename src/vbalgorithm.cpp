#include <RcppArmadillo.h>
#include "updateExpected.h"
#include "updateElbo.h"
#include "aux_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]
void vbalgorithm( const arma::mat& X, int D, int I, int J, int maxIter, double tolerance,
                  bool updatetau, std::string priorvar, arma::vec alphatau,
				          arma::vec betatau, arma::vec gammatau, arma::vec deltatau,
				          bool SVS, arma::vec priorInclusion, arma::vec beta1pi,
				          arma::vec beta2pi, double v0, std::string hypertype,
				          bool commonpi, int JD, arma::mat Tau, double qz,
				          bool scaleprior, Rcpp::List hpdis, int &it, double sigma2,
				          double aPostSigma, double denomX, double normX, double elbo,
				          double hW, double &finalElbo, double EWtauW, bool &converged,
				          arma::vec betastar1, arma::vec betastar2, arma::mat muW, arma::mat muP,
				          arma::mat W2, const arma::mat& XTX, arma::mat incProbs, arma::mat deltataupost,
                  double &globalElbo, arma::mat &globalMuW, arma::mat &globalMuP,
                  arma::mat &globaltau, double &globalSigma2,
                  arma::mat &globalbetatau, bool &globalconverged, arma::vec &elbovals,
                  Rcpp::List &globalhpdi, arma::vec &globalPriorInc,
                  arma::mat &globalIncPr, bool globalvar, bool hpdi, bool verbose ){

	finalElbo = 0.0;
	arma::vec allElbos;
  	arma::vec elboValue = arma::vec( maxIter + 1);


	elboValue.fill( -(arma::datum::inf) );
	converged = FALSE;

	for( it = 1; it < (maxIter + 1); it++ ){
		updateExpected( muW, W2, sigma2, muP, D, J, I, hW, X, XTX, Tau,
						hpdi, hpdis, qz, scaleprior, EWtauW, v0, incProbs, SVS,
						denomX, aPostSigma, normX,
						priorvar,
						updatetau, globalvar, gammatau, alphatau,
						betatau, hypertype, deltataupost, deltatau, JD,
						priorInclusion, betastar1, betastar2, commonpi,
						beta1pi, beta2pi );

		updateElbo( elbo, sigma2, globalvar, priorvar, W2, SVS, incProbs, v0,
					gammatau, betatau, deltataupost, alphatau, deltatau,
					hypertype, J, D, Tau, I, denomX,
					hW, JD, commonpi, beta1pi,
					beta2pi, priorInclusion, betastar1, betastar2, scaleprior );

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



		if( priorvar == "invgamma" ){

			if( all( gammatau > 0.0 ) ){

				if( hypertype == "common" ){

					globalbetatau.set_size(1, D);
					globalbetatau.row(0) = ((double(JD)*alphatau.t()) + gammatau.t()) / deltataupost(0,0);

				}else if( hypertype == "component" ){

					globalbetatau.set_size(1,D);
					globalbetatau.row(0) = ((double(J)*alphatau.t()) + gammatau.t()) / deltataupost.row(0);

				}else{

					globalbetatau.set_size(J,D);
					globalbetatau = 1 / deltataupost;
					globalbetatau.each_row() %= (alphatau.t() + gammatau.t());

				}

			}

		}

		globalPriorInc = priorInclusion;
		globalIncPr = incProbs;

	}

}
