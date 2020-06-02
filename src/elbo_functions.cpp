#include <RcppArmadillo.h>
#include "aux_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]
/***********************************************
*		Observed Data Log-Likelihood	      *
***********************************************/
void loglikelihood( double &loglik, double sigma2, int J, int I, double denomX ){
	loglik = (0.5 * (double(J) * double(I)) * ( -log(2.0 * PI * sigma2) ) ) - ( 0.5 * (1 / sigma2) * denomX );
}

/****************************
*		Prior for W 		 				*
****************************/
void logPriorW( double &logW, arma::mat Tau, int JD, arma::mat W2, int J,
					      int D, arma::mat logvar ){
	logW = arma::accu( (0.5 * logvar) - (0.5 * ( W2 % ( Tau ) ) )) - double(JD) * 0.5 * log( 2.0 * PI ) ;
}

/********************************
* 	Prior/Entropy for Tau      	*
********************************/
void priorentropyTau( double &logtau, double &hTau,
            		  int J, int D, arma::mat logvar, bool globalvar,
      				  arma::mat f, arma::mat Tau,
   					  arma::vec alphatau, arma::vec betatau, int JD ){
	int d;

	if( !globalvar ){ 															// Local prior variances
		for( int j=0; j<J; j++ ){
			for( d=0; d<D; d++ ){
				logtau += gamd( Tau(j,d), logvar(j,d), alphatau(d), betatau(d) );
				hTau += gamh( alphatau(d) + 0.5, exp( f( j,d ) ) ) ;
			}
		}
	}else{																		// Global prior variances
		for( d=0; d<D; d++ ){
			logtau += gamd( Tau(0,d), logvar(0,d), alphatau(d), betatau(d) );
			hTau += gamh( alphatau(d)+double(J)*0.5, exp( f( 0,d ) ) ) ;
		}
	}
}
