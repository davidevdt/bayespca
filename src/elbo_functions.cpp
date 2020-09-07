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
*		Prior for W 		 *
****************************/
void logPriorW( double &logW, arma::mat Tau, int JD, arma::mat W2, int J,
				int D, arma::mat logvar, double invsigma, bool SVS,
				double v0, arma::mat incProbs ){

	arma::mat tmpExp = W2 % ( Tau );
	tmpExp *= invsigma;

	if( SVS ){
		logvar += (1.0 - incProbs) * log(1 / v0);
		//logvar += ( incProbs * log(v0) ) ; 
		tmpExp = tmpExp % ( ((1.0 - incProbs) * ( 1 / v0 )) + incProbs);
		//tmpExp = tmpExp % ( (1.0 - incProbs) + (v0 + incProbs) );
	}

	logvar += log( invsigma );

	logW = arma::accu( (0.5 * logvar) - (0.5 * tmpExp)) - double(JD) * 0.5 * log( 2.0 * PI ) ;

}




/********************************
* 	Prior/Entropy for tau      	*
********************************/
void priorentropyTau( double &logtau, double &hTau, double &priorb, double &hb,
					  					int J, int D, arma::mat logvar, bool globalvar,
                      std::string priorvar, arma::mat f, arma::mat Tau,
                      arma::vec alphatau, arma::vec betatau, arma::vec gammatau,
                      arma::vec deltatau, arma::mat deltataupost, std::string hypertype, int JD ){

	int d;


	if( !globalvar ){ 															// Local prior variances
		if( priorvar == "jeffrey" ){
			// logtau = arma::accu( -logvar );
			logtau = arma::accu( logvar );
			for( int j=0; j<J; j++ ){
				for( d=0; d<D; d++ ){
					hTau += gamh( 0.5, exp( f( j,d ) ) ) ;
				}
			}
		}else{																	// priorvar == "invgamma"
			if( any(gammatau <= 0.0) ){
				for( int j=0; j<J; j++ ){
					for( d=0; d<D; d++ ){
						logtau += gamd( Tau(j,d), logvar(j,d), alphatau(d), betatau(d) );
						hTau += gamh( alphatau(d) + 0.5, exp( f( j,d ) ) ) ;
					}
				}
			}else{
				if( hypertype == "common" ){
					double btmp = (gammatau(0) + (double(JD)*alphatau(0))) / deltataupost(0,0);
					priorb += gamd( btmp, (diGammaFunc( (double(JD)*alphatau(0)) + gammatau(0) )) - log( deltataupost(0,0) ), gammatau(0), deltatau(0) );
					hb += gamh( (double(JD)*alphatau(0)) + gammatau(0), deltataupost(0,0) );
					for( d=0; d<D; d++ ){
						for( int j=0; j<J; j++ ){
							logtau += gamd( Tau(j,d), logvar(j,d), alphatau(0), btmp );
							hTau += gamh( alphatau(0)+0.5, exp( f( j,d ) ) ) ;
						}
					}
				}else if( hypertype == "component" ){
					double btmp;
					for( d=0; d<D; d++ ){
						btmp = (gammatau(d) + (double(J)*alphatau(d))) / deltataupost(0,d);
						priorb += gamd( btmp, (diGammaFunc( (double(J)*alphatau(d)) + gammatau(d) )) - log( deltataupost(0,d) ), gammatau(d), deltatau(d) );
						hb += gamh( (double(J)*alphatau(d)) + gammatau(d), deltataupost(0,d) );
						for( int j=0; j<J; j++ ){
							logtau += gamd( Tau(j,d), logvar(j,d), alphatau(d), btmp );
							hTau += gamh( alphatau(d)+0.5, exp( f( j,d ) ) ) ;
						}
					}
				}else{
					double btmp;
					for( int j=0; j<J; j++ ){
						for( d=0; d<D; d++ ){
							btmp = (gammatau(d) + alphatau(d)) / deltataupost(j,d);
							priorb += gamd( btmp, (diGammaFunc( alphatau(d) + gammatau(d) )) - log( deltataupost(j,d) ), gammatau(d), deltatau(d) );
							hb += gamh( alphatau(d) + gammatau(d), deltataupost(j,d) );
							logtau += gamd( Tau(j,d), logvar(j,d), alphatau(d), btmp );
							hTau += gamh( alphatau(d)+0.5, exp( f( j,d ) ) );
						}
					}
				}
			}
		}
	}else{																		// Global prior variances

		if( priorvar == "jeffrey" ){
			logtau = arma::accu( -logvar.row(0) );
			for( d=0; d<D; d++ ){
				hTau += invgamh( double(J)*0.5, exp( f( 0,d ) ) );
			}
		}else{																	// priorvar == "invgamma"
			if( any( gammatau <= 0.0 ) ){									    // Fixed Hyperparameters
				for( d=0; d<D; d++ ){
					logtau += gamd( Tau(0,d), logvar(0,d), alphatau(d), betatau(d) );
					hTau += gamh( alphatau(d)+double(J)*0.5, exp( f( 0,d ) ) ) ;
				}
			}else{ 										// Random ( component-specific ) hyperparameters
				double btmp;
				for( d=0; d<D; d++ ){
					btmp = (gammatau(d) + (alphatau(d))) / deltataupost(0,d);
					priorb += gamd( btmp, (diGammaFunc( (alphatau(d)) + gammatau(d) ))-log( deltataupost(0,d) ), gammatau(d), deltatau(d) );
					hb += gamh( (alphatau(d)) + gammatau(d), deltataupost(0,d) );
					for( int j=0; j<J; j++ ){
						logtau += gamd( Tau(0,d), logvar(0,d), alphatau(d), btmp );
						hTau += gamh( alphatau(d)+(double(J)*0.5), exp( f( 0,d ) ) ) ;
					}
				}
			}
		}
	}
}


/********************************************
*	Prior/Entropy for SVS components	   *
********************************************/
void priorentropySVS( double &logPriorIncProbs, double &hPriorIncProbs, double &logPriorGlobalProb,
					  double &hPriorGlobalProb, bool commonpi, arma::vec beta1pi,
					  arma::vec beta2pi, arma::vec priorInclusion, arma::vec betastar1,
					  arma::vec betastar2, arma::mat incProbs, int D ){

	arma::vec logpg(D);
	arma::vec logpinvg(D);

	if( commonpi ){
		if( any(beta1pi <= 0.0 )){
			logpg.fill( log( priorInclusion(0) ) );
			logpinvg.fill( log( (1.0)-priorInclusion(0) ) );
		}else{
			double sumtmp = betastar1(0) + betastar2(0);
			logpg.fill( diGammaFunc(betastar1(0)) - diGammaFunc(sumtmp) );
			logpinvg.fill( diGammaFunc(betastar2(0)) - diGammaFunc(sumtmp) );
		}

		logPriorIncProbs += logpg(0) * arma::accu( incProbs ) + logpinvg(0) * arma::accu( 1.0 - incProbs );

	}else{
		int d;
		if( any(beta1pi <= 0.0) ){
			for( d=0; d<D; d++ ){
				logpg(d) = log(priorInclusion(d));
				logpinvg(d) = log( (1.0)-priorInclusion(d) );
				logPriorIncProbs += logpg(d) * arma::accu( incProbs.col(d) ) + logpinvg(d) * arma::accu( 1.0 - incProbs.col(d) );
			}
		}else{
			double sumtmp;
			for( d=0; d<D; d++ ){
				sumtmp = betastar1(d) + betastar2(d);
				logpg(d) = diGammaFunc(betastar1(d)) - diGammaFunc(sumtmp);
				logpinvg(d) = diGammaFunc(betastar2(d)) - diGammaFunc(sumtmp);
				logPriorIncProbs += logpg(d) * arma::accu( incProbs.col(d) ) + logpinvg(d) * arma::accu( 1.0 - incProbs.col(d) );
			}
		}
	}

	arma::mat hPrIncl = ( (1.0 - incProbs) % log(1.0 - incProbs) + (incProbs) % log(incProbs) );
	hPrIncl.elem( find_nonfinite(hPrIncl) ).zeros();
	hPriorIncProbs -= arma::accu( hPrIncl );

	if( all(beta1pi > 0.0) ){
		if( commonpi ){
			logPriorGlobalProb += (beta1pi(0)-1.0)*(logpg(0)) + (beta2pi(0)-1.0)*(logpinvg(0)) - betaFunc( beta1pi(0), beta2pi(0), TRUE );
			hPriorIncProbs += betaFunc( betastar1(0), betastar2(0), TRUE ) - (betastar1(0)-1.0)*diGammaFunc( betastar1(0) ) - (betastar2(0)-1.0)*diGammaFunc( betastar2(0) ) + (betastar1(0)+betastar2(0)-2.0)*diGammaFunc( betastar1(0) + betastar2(0) );
		}else{
			int d;
			for( d=0; d<D; d++ ){
				logPriorGlobalProb += (beta1pi(d)-1.0)*(logpg(d)) + (beta2pi(d)-1.0)*(logpinvg(d)) - betaFunc( beta1pi(d), beta2pi(d), TRUE );
				hPriorIncProbs += betaFunc( betastar1(d), betastar2(d), TRUE ) - (betastar1(d)-1.0)*diGammaFunc( betastar1(d) ) - (betastar2(d)-1.0)*diGammaFunc( betastar2(d) ) + (betastar1(d)+betastar2(d)-2.0)*diGammaFunc( betastar1(d) + betastar2(d) );
			}
		}
	}
}
