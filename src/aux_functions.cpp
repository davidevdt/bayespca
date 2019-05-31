#include <RcppArmadillo.h>
#include <math.h>         
#include <float.h>

// [[Rcpp::depends(RcppArmadillo)]]

#define PI M_PI 


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
arma::mat invSMW( arma::mat D, arma::mat X, int I ){
	return D - D * X.t() * arma::inv( arma::eye(I, I) + X * D * X.t() ) * X * D;
}



/********************************
*	Extract trace of X	 	*
********************************/
void funcX( arma::mat X, arma::mat &XTX,  double &trXTX ){
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
void updateIncProbs( arma::mat &incProbs, int J, int D, arma::vec logpg, 
                     arma::vec logpinvg, arma::mat W2, double v0, 
                     arma::mat invTau, double invSigma2 ){
	int j; 
	int d; 
	
	arma::vec valtmp(J); 
	
	double ptmp; 
	double odds; 
	
	for( d=0; d<D; d++ ){
		odds = (logpg(d) - logpinvg(d)) + (0.5 * log( v0 ));  
 
		valtmp = (W2.col(d) * invSigma2) % ( invTau.col(d) );
		for( j=0; j<J; j++ ){				
			ptmp =  - 0.5 * (valtmp(j)) * (1.0 - (1/v0));
			incProbs(j,d) = 1.0 / (1.0 + exp( -(ptmp + odds) ) );
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
				f = f % (((1.0-incProbs)*(1/v0)) + incProbs); 
			}
      
			f = arma::log(f);
		}else{ 																// priorvar == "invgamma" 
			f = W2 * invsigma / 2.0; 
      
			if( SVS ){
				f = f % (((1.0-incProbs)*(1/v0)) + incProbs); 
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
					 arma::mat invTau, std::string hypertype, arma::vec alphatau 
                    ){
  
	arma::mat logvar(J,D); 
  
	if( !globalvar ){
    
		if( priorvar == "fixed" ){
			logvar = log( 1 / invTau ); 
		}else if( priorvar == "jeffrey" ){
			logvar = f - diGammaFunc( 0.5 ) ; 
		}else{
			int d;
			logvar = f; 
			for( d=0; d<D; d++ ){
			logvar.col(d) -= diGammaFunc( alphatau(d) + 0.5 );
			}
		}
	
	}else{     																 // globalvar = TRUE
		if( priorvar == "fixed" ){
			logvar = arma::log( 1/invTau );  
		}else if( priorvar == "jeffrey" ){
			logvar.each_row() = f.row(0) - diGammaFunc( double(J)*0.5 ) ;  
		}else{  																// priorvar == "invgamma"
			int d;
			logvar = f; 
			for( d=0; d<D; d++ ){
				logvar.col(d) -= diGammaFunc( alphatau(d) + (double(J)*0.5));
			}      
		}
	}	
  
	logvar.elem(find_nonfinite(logvar)).zeros();
  
	return logvar;   
}


