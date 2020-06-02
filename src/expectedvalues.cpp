#include <RcppArmadillo.h>
#include "aux_functions.h"

#define MAX_FLOAT 1E+37

// [[Rcpp::depends(RcppArmadillo)]]

/****************************
*	 Update W, P, sigma2 	*
****************************/
void updateWPXsigma( arma::mat &muW, arma::mat &W2, double &sigma2, arma::mat &muP, int D,
                     int J, int I, double &hW, const arma::mat& X, const arma::mat& XTX, arma::mat Tau,
                     std::string priorvar, bool hpdi, Rcpp::List &hpdis, double qz,
                     double &denomX, double aPostSigma, double normX ){

	int d;

	double Det;
	double tmpDet;
	double invSigma2 = 1.0 / sigma2;

	arma::vec taucol(J);

	arma::mat sigmaWd(J,J);
	arma::mat id = arma::eye(J,J);
	arma::mat EWWT(J,J);

	hW = 0.0;
	EWWT.zeros();

  // 1. Update W
	muW.zeros();

	for( d = 0; d<D; d++ ){
		taucol = Tau.col( d );
	  sigmaWd = arma::inv_sympd(arma::diagmat(taucol) + (invSigma2 * XTX ));

		muW.col( d ) = invSigma2 * (sigmaWd * XTX * muP.col(d) );
		W2.col( d ) = ( arma::square( muW.col( d ) ) + sigmaWd.diag() );

		// Entropy for W
		tmpDet = arma::det(sigmaWd);
		Det = (tmpDet > MAX_FLOAT) ? MAX_FLOAT : tmpDet ;
		if( Det != 0.0 ){
		  hW += 0.5 * log( ( Det ) ) + ((double(J) / 2.0) * ( 1.0 + log( 2.0 * PI ) )) ;
		}else{
			hW += 0.0;
		}

		EWWT += sigmaWd;

		if( hpdi ){
			hpdis( d ) = retHPDI( muW.col(d), arma::sqrt(sigmaWd.diag()), qz, J );
		}
	}

	EWWT += muW * muW.t();


	// 2. Update P
	arma::mat U(J, D);
	arma::vec DG(D);
	arma::mat V(J, D);
	SVD( XTX * muW, U, DG, V, D, D );
	muP = U * V.t();

	// 3. Update sigma2
	denomX = normX - 2 * arma::trace( muW.t() * XTX * muP ) + arma::trace( XTX * EWWT );
	double bPostSigma = denomX;
	sigma2 =  bPostSigma / aPostSigma;

}


/******************************
*         Update Tau          *
 *****************************/
 void updateTau( std::string priorvar, bool updatetau,
                 int D, arma::mat &Tau, arma::mat W2, bool globalvar,
                 arma::vec alphatau, arma::vec betatau, int JD, int J ){

	int d;

	if( !globalvar ){      							   // Local prior variances

	 if( priorvar == "fixed" ){        		// Fixed Tau
				for( d=0; d<D; d++ ){
					Tau.col(d) =  1. / ( W2.col(d) ) ;
				}
		}else{      		                    // Gamma priors
				for( d=0; d<D; d++ ){
					Tau.col(d) = (alphatau(d) + 0.5) / ( betatau(d) + ( W2.col(d) / (2.0) ) )  ;
				}
		}

  }else{															 // Global prior variances

		if( priorvar == "fixed" ){         // Fixed tau
				for( d=0; d<D; d++ ){
					Tau.col(d).fill( double(J) / arma::accu(W2.col(d)) ) ;
				}
		}else{                              // Gamma priors
			for( d=0; d<D; d++ ){
				Tau.col(d).fill( (alphatau(d) + ( double(J) * 0.5 ) ) / ( betatau(d) + arma::accu( W2.col(d) / (2.0) ) )) ;
			}
		}
	}

}
