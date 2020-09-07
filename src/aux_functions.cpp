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
*		Gamma, logGamma, Digamma Functions, Beta	  *
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



/********************
*		SVS		  *
********************/
void updateIncProbs( arma::mat &incProbs, int J, int D, arma::vec priorInclusion,
                     arma::mat W2, double v0,
                     arma::mat Tau, double invSigma2 ){
	int j;
	int d;

	double p_inc_prop;
	double p_exc_prop;

	for( d=0; d<D; d++ ){
		for( j=0; j<J; j++ ){
			p_inc_prop = priorInclusion(d) * exp( - 0.5 * Tau(j, d) * W2(j,d) * invSigma2 );
			p_exc_prop = ((1.0 - priorInclusion(d)) / sqrt(v0) ) * exp( - 0.5 * Tau(j, d) * W2(j,d) * (1.0 / v0) * invSigma2 );		
			incProbs(j,d) = p_inc_prop / (p_inc_prop + p_exc_prop);
		}
	}
}

/************************
*   FUNCTIONS FOR ELBO	 *
************************/
arma::mat fMat( bool globalvar, std::string priorvar, arma::mat W2, double invsigma, bool SVS,
                arma::mat incProbs, double v0, arma::vec gammatau, arma::vec betatau,
                arma::mat deltataupost, std::string hypertype,
                int J, int D, arma::vec alphatau, int JD ){

	arma::mat f;

	if( !globalvar ){
		if( priorvar == "fixed" ){
			f.set_size(1,1);
		}else if( priorvar == "jeffrey" ){
			f = W2 * (invsigma / 2.0) ;

			if( SVS ){
				// f = f % (((1.0-incProbs)*(1/v0)) + incProbs);
				f = f % ((1.0-incProbs) + (v0 * incProbs));
			}

			f = arma::log(f);
		}else{ 																// priorvar == "invgamma"
			f = W2 * invsigma / 2.0;

			if( SVS ){
				f = f % (((1.0-incProbs)*(1/v0)) + incProbs);
				// f = f % ((1.0-incProbs) + (v0 * incProbs));
			}

			if( any( gammatau <= 0.0 ) ){
				f.each_row() += betatau.t();
				f = log( f );
			}else{
				if( hypertype == "common" ){
					double btmp = ( gammatau(0) + (double(JD)*alphatau(0)) ) / deltataupost(0,0);
					f += btmp;
					f = log( f );
				}else if( hypertype == "component" ){
					arma::vec btmp = ( gammatau + (double(J)*alphatau) ) / deltataupost.row(0).t();
					f.each_row() += btmp.t();
					f = log( f );
				}else{
					arma::mat btmp = 1 / deltataupost;
					btmp.each_row() %= ( gammatau.t() + alphatau.t() );
					f += btmp;
					f = log(f);
				}
			}
		}
	}else{    																	// globalvar
		if( priorvar == "fixed" ){
			f.set_size(1,1);
		}else if( priorvar == "jeffrey" ){
			f.set_size(1,D);
			if( !SVS ){
				f.each_row() = arma::sum( W2, 0 );
				f *= invsigma / 2.0;
				f = arma::log(f);
			}else{
				f.each_row() = arma::sum( W2 % ( ((1.0-incProbs)*(1.0/v0)) + incProbs), 0 );
				// f.each_row() = arma::sum( W2 % ((1.0-incProbs) + (v0 * incProbs)), 0 );
				f *= invsigma / 2.0;
				f = arma::log(f);
			}
		}else{    																// priorvar == "invgamma"
			f.set_size(J,D);
			if( any( gammatau <= 0.0 ) ){
				if( !SVS ){
					f.each_row() = arma::sum( W2, 0 );
					f *= invsigma / 2.0;
				}else{
					f.each_row() = arma::sum( W2 % ( ((1.0-incProbs)*(1.0/v0)) + incProbs ), 0 );
					// f.each_row() = arma::sum( W2 % ((1.0-incProbs) + (v0 * incProbs)), 0 );
					f *= invsigma / 2.0;
				}
				f.each_row() += betatau.t();
				f = arma::log(f);
			}else{
				arma::vec btmp = (gammatau + (alphatau)) / deltataupost.row(0).t();
				if( !SVS ){
					f.each_row() = arma::sum( W2, 0 );
					f *= invsigma / 2.0;
				}else{
					f.each_row() = arma::sum( W2 % ( ((1.0-incProbs)*(1.0/v0)) + incProbs ), 0 );
					// 	f.each_row() = arma::sum( W2 % ((1.0-incProbs) + (v0 * incProbs)), 0 );
					f *= invsigma / 2.0;
				}
				f.each_row() += btmp.t();
				f = arma::log( f );
			}
		}
	}
	return f;
}


arma::mat logvarMat( bool globalvar, int J, int D, arma::mat f, std::string priorvar,
					 arma::mat Tau, std::string hypertype, arma::vec alphatau
                    ){

	arma::mat logvar(J,D);

	if( !globalvar ){

		if( priorvar == "fixed" ){
			logvar = log( Tau );
		}else if( priorvar == "jeffrey" ){
			logvar = diGammaFunc( 0.5 ) - f ;
		}else{
			int d;
			logvar = - f;
			for( d=0; d<D; d++ ){
			logvar.col(d) += diGammaFunc( alphatau(d) + 0.5 );
			}
		}

	}else{     																 // globalvar = TRUE
		if( priorvar == "fixed" ){
			logvar = arma::log( Tau );
		}else if( priorvar == "jeffrey" ){
			logvar.each_row() = diGammaFunc( double(J)*0.5 ) - f.row(0);
		}else{  																// priorvar == "invgamma"
			int d;
			logvar = -f;
			for( d=0; d<D; d++ ){
				logvar.col(d) += diGammaFunc( alphatau(d) + (double(J)*0.5));
			}
		}
	}

	logvar.elem(find_nonfinite(logvar)).zeros();

	return logvar;
}
