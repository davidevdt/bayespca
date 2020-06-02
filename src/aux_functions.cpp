#include <RcppArmadillo.h>
#include <math.h>
#include <float.h>

// [[Rcpp::depends(RcppArmadillo)]]

#define PI M_PI

/********************************************************
* 		Compute SVD Decomposition of a matrix M 		******
********************************************************/
Rcpp::Function svdExt("svd");

void SVD( const arma::mat&M, arma::mat& U, arma::vec& D, arma::mat& V, int nu, int nv ){
	Rcpp::List svdList = svdExt(M, Rcpp::Named("nu") = nu, Rcpp::Named("nv") = nv, Rcpp::Named("LINPACK") = false );

	U = Rcpp::as<arma::mat>(svdList["u"]);
	D = Rcpp::as<arma::vec>(svdList["d"]);
	V = Rcpp::as<arma::mat>(svdList["v"]);
}

/******************************************************
*		Gamma, logGamma, Digamma, Beta Functions		  		*
******************************************************/
double gammaFunc( double x, bool logScale ){
	double ret;
	if( !logScale ){
		if( x == 0.5 ){
			ret = sqrt( PI );
		}else{
			ret = R::gammafn( x );
		}
	}else{
		ret = R::lgammafn( x );
	}
	return ret;
}

double diGammaFunc( double x ){
	double ret;
	ret = R::digamma( x );
	return ret;
}

double betaFunc( double a, double b, bool logScale ){
	double ret;
	ret = gammaFunc( a, TRUE ) + gammaFunc( b, TRUE ) - gammaFunc( (a+b), TRUE );
	if( logScale ){
		return ret;
	}else{
		return exp( ret );
	}
}

/**************************************************
*	Sherman - Morrison - Woodbury Inversion 		*
**************************************************/
// Note: this function is not used in new package versions.

// arma::mat invSMW( arma::mat D, arma::mat X, int I ){
//	return D - D * X.t() * arma::inv( arma::eye(I, I) + X * D * X.t() ) * X * D;
// }

/********************************
*	Extract trace of X	 	*
********************************/
void funcX( const arma::mat& X, arma::mat &XTX,  double &trXTX ){
	XTX = X.t() * X;
	trXTX = arma::trace( XTX );
}

/****************
*	Densities  *
****************/
// log-density of Gamma distribution
double gamd( double x, double logx, double a, double b ){
	return a * log(b) + (a - 1.0) * logx - b * x - gammaFunc( a, TRUE );
}

// log-density of Inverse-Gamma distribution
double invgamd( double invx, double logx, double a, double b ){
	return a * log(b) - (a + 1.0) * logx - b * invx - gammaFunc( a, TRUE );
}

// Entropy of Gamma distribution
double gamh( double a, double b ){
	return a - log( b ) + gammaFunc( a, TRUE ) + (1.0 - a) * diGammaFunc(a);
}

// Entropy of Inverse Gamma distribution
double invgamh( double a, double b ){
	return a + log( b * exp( gammaFunc( a, TRUE ) ) ) - (1.0 + a) * diGammaFunc(a);
}

/************************
*	HPDIs calculation	*
************************/
arma::mat retHPDI( arma::vec mu, arma::vec sigma, double qz, int J ){
	arma::mat ret(J, 2, arma::fill::zeros);
	int j;

	for( j = 0; j < J; j++ ){
		ret( j, 0 ) = mu(j) - (qz * sigma(j));
		ret( j, 1 ) = mu(j) + (qz * sigma(j));
	}

	return ret;
}

/************************
*   FUNCTIONS FOR ELBO	 *
************************/
arma::mat fMat( bool globalvar, std::string priorvar, arma::mat W2,
                arma::vec betatau, int J, int D, arma::vec alphatau, int JD ){

	arma::mat f;

	if( !globalvar ){
		if( priorvar == "fixed" ){
			f.set_size(1,1);
		}else{ 																// priorvar == "gamma"
			f = W2 / 2.0;
			f.each_row() += betatau.t();
			f = log( f );
		}
	}else{    																	// globalvar
		if( priorvar == "fixed" ){
			f.set_size(1,1);
		}else{    																// priorvar == "gamma"
			f.set_size(J,D);
			f.each_row() = arma::sum( W2, 0 );
			f = f / 2.0;
			f.each_row() += betatau.t();
			f = arma::log(f);
		}
	}
	return f;
}


arma::mat logvarMat( bool globalvar, int J, int D, arma::mat f, std::string priorvar,
					 					 arma::mat Tau, arma::vec alphatau ){

	arma::mat logvar(J,D);

	if( !globalvar ){
		if( priorvar == "fixed" ){
			logvar = log( Tau );
		}else{																	// priorvar == "gamma"
			int d;
			logvar = - f;
			for( d=0; d<D; d++ ){
			logvar.col(d) += diGammaFunc( alphatau(d) + 0.5 );
			}
		}
	}else{     																 // globalvar
		if( priorvar == "fixed" ){
			logvar = arma::log( Tau );
		}else{  																// priorvar == "gamma"
			int d;
			logvar = -f;
			for( d=0; d<D; d++ ){
				logvar.col(d) += diGammaFunc( alphatau(d) + (double(J) * 0.5));
			}
		}
	}

	logvar.elem(find_nonfinite(logvar)).zeros();
	return logvar;
}
