#include <RcppArmadillo.h>
#include "aux_functions.h"

#define MAX_FLOAT 1E+37

// [[Rcpp::depends(RcppArmadillo)]]


/****************************  
***	Update W, P, sigma2 	***
****************************/
void updateWPXsigma( arma::mat &muW, arma::mat &W2, double &sigma2, arma::mat &muP, int D, 
                     int J, int I, double &hW, const arma::mat& X, const arma::mat& XTX, arma::mat invTau, 
                     std::string priorvar, bool hpdi, Rcpp::List &hpdis, double qz, 
                     bool scaleprior, double &EWtauW, double v0, 
                     arma::mat incProbs, double SVS, double &denomX, 
                     double aPostSigma, double normX ){
	
	
	int d; 
	double Det; 
	double tmpDet; 
	double invSigma2 = 1.0 / sigma2; 
	double invsigma = sqrt(1 / sigma2);

	arma::mat sigmaWd(J,J); 
	arma::mat id = arma::eye(J,J);
	arma::mat taucol(J,J);		
	arma::mat EWWT(J,J); 
	
	hW = 0.0; 
	
	if( scaleprior ){
		EWtauW = 0.0;
	}
	
	// Update W 
	EWWT.zeros(); 
	muW.zeros();
	
	for( d = 0; d<D; d++ ){		

		taucol = arma::inv( arma::diagmat( invTau.col( d ) ) );

		if( scaleprior ){
			taucol *= (sigma2); 
		}
		
		if( SVS ){
			taucol = taucol % arma::diagmat( v0 / ( (1.0 - incProbs.col(d)) + ((incProbs.col(d)) * v0) ) );  
		}
	 
		if( double(J / I) <= 100. ){
			sigmaWd = arma::inv(id +  (taucol  * (invSigma2 * XTX ))) * taucol ;
		}else{
			sigmaWd = invSMW( taucol, invsigma * X, I );
		}
	
		muW.col( d ) = invSigma2 * ( sigmaWd * XTX * muP.col(d) ); 
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

		if( scaleprior ){
			taucol = arma::diagmat( invTau.col( d ) );
			if( SVS ){
				taucol = taucol % arma::diagmat( ( (1.0 - incProbs.col(d) ) * (1.0 / v0) ) + incProbs.col(d) );
			}
			EWtauW += arma::accu( muW.col(d).t() * taucol * muW.col(d) ) + arma::trace( arma::sqrt(taucol) * sigmaWd * arma::sqrt(taucol.t()) );
		} 

	}	
	
	EWWT += muW * muW.t(); 

	//	Update P 
	arma::mat U(J, D); 
	arma::vec DG(D); 
	arma::mat V(J, D); 	

	SVD( XTX * muW, U, DG, V, D, D );			
	muP = U * V.t();	

	// Update sigma2 
	denomX = normX - 2 * arma::trace( muW.t() * XTX * muP ) + arma::trace( XTX * EWWT ); 
	double bPostSigma = denomX; 

	if( scaleprior ){
		bPostSigma += ( EWtauW ); 
	}
	
	sigma2 =  bPostSigma / aPostSigma; 

}


/******************************
*         Update tau          * 
 *****************************/
void updateTau( double sigma2, bool scaleprior, std::string priorvar, bool updatetau,
                   bool SVS, int D, arma::mat &invTau, arma::mat W2, bool globalvar,
                   arma::mat incProbs, double v0, arma::vec gammatau, arma::vec alphatau,
                   arma::vec &betatau, std::string hypertype, arma::mat &deltataupost,
                   arma::vec deltatau, int JD, int J
               ){
  
	int d;
	double tmpinvsigma; 
  
	if( scaleprior ){
		tmpinvsigma = 1.0/sigma2; 
	}else{
		tmpinvsigma = 1.0; 
	}
  
	if( !globalvar ){      										// Local prior variances
		
		// Fixed tau
		if( priorvar == "fixed" ){
			if( !SVS ){
				for( d=0; d<D; d++ ){			
					invTau.col(d) =  ( 1.0 / (W2.col(d) * tmpinvsigma) ) ; 
				}
			}else{
				for( d=0; d<D; d++ ){
					invTau.col(d) =  ( 1.0 / ((W2.col(d) * tmpinvsigma) % ( ((1.0-incProbs.col(d))*(1.0/v0)) + incProbs.col(d) ) ) ); 
				}
			}
		// Jeffrey's prior 
		}else if( priorvar == "jeffrey" ){
			if( !SVS ){
				for( d=0; d<D; d++ ){
					invTau.col(d) = 1.0 / (W2.col(d) * tmpinvsigma ); 
				}
			}else{			
				for( d=0; d<D; d++ ){
					invTau.col(d) =  ( 1.0 / ((W2.col(d) * tmpinvsigma) % ( ((1.0-incProbs.col(d))*(1.0/v0)) + incProbs.col(d) ) ) );  
				}
			}
		// InverseGamma prior 
		}else{
			// Fixed hyperparameters 
			if( any(gammatau <= 0.0) ){
				if( !SVS ){
					for( d=0; d<D; d++ ){
						invTau.col(d) = (alphatau(d) + 0.5) / ( betatau(d) + ( W2.col(d) * tmpinvsigma / (2.0) ) ) ; 
					}
				}else{
					for( d=0; d<D; d++ ){
						invTau.col(d) = (alphatau(d) + 0.5) / ( betatau(d) + ( ((W2.col(d) * tmpinvsigma) % ( ((1.0-incProbs.col(d))*(1.0/v0)) + incProbs.col(d)) )/ (2.0) ) ) ; 
					}		
				}
			// Random beta
			}else{
				if( hypertype == "common" ){
					deltataupost(0,0) = deltatau(0) + arma::accu( invTau );	
					double btmp = (gammatau(0) + (double(JD)*alphatau(0))) / deltataupost(0,0); 
					if( !SVS ){
						for( d=0; d<D; d++ ){
							invTau.col(d) = (alphatau(0) + 0.5) / (btmp + ((W2.col(d) * tmpinvsigma ) / 2.0) ); 
						}
					}else{
						for(d=0; d<D; d++){
							invTau.col(d) = (alphatau(0) + 0.5) / (btmp + ( ((W2.col(d) * tmpinvsigma ) % ( ((1.0-incProbs.col(d))*(1.0/v0)) + incProbs.col(d)) )/ 2.0) );           
						}
					}
				}else if( hypertype == "component" ){
					deltataupost.each_row() = ( deltatau.t() + arma::sum( invTau, 0) );  
					double btmp;
					if( !SVS ){
						for( d=0; d<D; d++ ){
							btmp = (gammatau(d) + (double(J)*alphatau(d))) / deltataupost(0,d);
							invTau.col(d) = (alphatau(d) + 0.5) / (btmp + ((W2.col(d) * tmpinvsigma ) / 2.0) );
						}
					}else{
						for( d=0; d<D; d++ ){
							btmp = ( gammatau(d) + (double(J)*alphatau(d)) ) / deltataupost(0,d);
							invTau.col(d) = (alphatau(d) + 0.5) / (btmp + ( ((W2.col(d) * tmpinvsigma ) % ( ((1.0-incProbs.col(d))*(1.0/v0)) + incProbs.col(d)) )/ 2.0) );                       
						}
					}
				}else{													// hypertype == "local"
					int j; 
					double btmp; 
					if( !SVS ){
						for( j=0;j<J;j++ ){
							for( d=0; d<D; d++){
								deltataupost(j,d) = deltatau(d) + invTau(j,d);
								btmp = (gammatau(d) + alphatau(d)) / deltataupost(j,d);
								invTau(j,d) = (alphatau(d) + 0.5) / (btmp + ((W2(j,d)*tmpinvsigma)/2.0));
							}
						}
					}else{
						for( j=0; j<J; j++ ){
							for( d=0; d<D; d++ ){
								deltataupost(j,d) = deltatau(d) + invTau(j,d);
								btmp = (gammatau(d) + alphatau(d)) / deltataupost(j,d);
								invTau(j,d) = (alphatau(d) + 0.5) / (btmp + ((W2(j,d)*tmpinvsigma*( ((1.0-incProbs(j,d))*(1.0/v0)) + incProbs(j,d)) )/2.0));
							}
						}
					}
				}
			}
		}
	}else{															// Global prior variances
		// Fixed tau  
		if( priorvar == "fixed" ){
			if( !SVS ){
				for( d=0; d<D; d++ ){			
					invTau.col(d).fill( double(J) / arma::accu(W2.col(d) * tmpinvsigma) ) ; 
				}
			}else{
				for( d=0; d<D; d++ ){
					invTau.col(d).fill( double(J) / arma::accu((W2.col(d) * tmpinvsigma) % ( ((1.0-incProbs.col(d))*(1.0/v0)) + incProbs.col(d) ) ) ); 
				}
			}   
		// Jeffrey's prior 
		}else if( priorvar == "jeffrey" ){
			if( !SVS ){
				for( d=0; d<D; d++ ){
					invTau.col(d).fill(double(J) / arma::accu(W2.col(d) * tmpinvsigma )); 
				}
			}else{			
				for( d=0; d<D; d++ ){
					invTau.col(d).fill( double(J) / arma::accu((W2.col(d) * tmpinvsigma) % ( ((1.0-incProbs.col(d))*(1.0/v0)) + incProbs.col(d) ) ) );  
				}
			} 
		// InverseGamma priors 
		}else{
			// Fixed hyperparamters 
			if( any( gammatau <= 0.0 ) ){   												
				if( !SVS ){
					for( d=0; d<D; d++ ){
						invTau.col(d).fill( (alphatau(d) + ( double(J) * 0.5 ) ) / ( betatau(d) + arma::accu( W2.col(d) * tmpinvsigma / (2.0) ) )) ; 
					}
				}else{
					for( d=0; d<D; d++ ){
						invTau.col(d).fill( (alphatau(d) + ( double(J) * 0.5 ) ) / ( betatau(d) + arma::accu( ((W2.col(d) * tmpinvsigma) % ( ((1.0-incProbs.col(d))*(1.0/v0)) + incProbs.col(d) ) )/ (2.0) ) )) ; 
					}		
				}
			// Random ( component-specific ) beta  	
			}else{   																			 
				deltataupost.each_row() = ( deltatau.t() +  invTau.row(0)  );
				double btmp;
				if( !SVS ){
					for( d=0; d<D; d++ ){
						btmp = (gammatau(d) + (alphatau(d))) / deltataupost(0,d);
						invTau.col(d).fill((alphatau(d) + (double(J)*0.5)) / (btmp + arma::accu((W2.col(d) * tmpinvsigma ) / 2.0) ));
					}
				}else{
					for( d=0; d<D; d++ ){
						btmp = (gammatau(d) + (alphatau(d))) / deltataupost(0,d);
						invTau.col(d).fill((alphatau(d) + (double(J)*0.5)) / (btmp + arma::accu(((W2.col(d) * tmpinvsigma ) % ( ((1.0-incProbs.col(d))*(1.0/v0)) + incProbs.col(d)) )/ 2.0) ) );                       
					}
				}        
			}
		}
	}
}



/******************
*	Update SVS	*
******************/
void updateSVS( arma::mat &incProbs, arma::vec &priorInclusion, arma::vec &betastar1, 
                   arma::vec &betastar2, int J, int D, arma::mat W2, double v0, 
                   arma::mat invTau, double sigma2, bool scaleprior, 
                   bool commonpi, arma::vec beta1pi, arma::vec beta2pi, int JD ){
	
	arma::vec logpg(D);
	arma::vec logpinvg(D); 
	
	
	double invSigma2;
	if( scaleprior  ){
		invSigma2 = 1/sigma2;
	}else{
		invSigma2 = 1.0; 
	} 
	
	if( commonpi ){
		if( any(beta1pi <= 0.0) ){
			logpg.fill( log( priorInclusion(0) ) );
			logpinvg.fill( log( (1.0)-priorInclusion(0) ) );
		}else{
			double sumtmp = betastar1(0) + betastar2(0);
			logpg.fill( diGammaFunc(betastar1(0)) - diGammaFunc(sumtmp) );
			logpinvg.fill( diGammaFunc(betastar2(0)) - diGammaFunc(sumtmp) );
		}
	}else{
		int d; 
		if( any( beta1pi <= 0.0) ){
			for( d=0; d<D; d++ ){
				logpg(d) = log( priorInclusion(d) );
				logpinvg(d) = log( (1.0)-priorInclusion(d) );			
			}
		}else{
			double sumtmp; 
			for( d=0; d<D; d++ ){
				sumtmp = betastar1(d) + betastar2(d); 		
				logpg(d) = diGammaFunc(betastar1(d)) - diGammaFunc(sumtmp);
				logpinvg(d) = diGammaFunc(betastar2(d)) - diGammaFunc(sumtmp);				
			}			
		}		
	}
	
	updateIncProbs( incProbs, J, D, logpg, logpinvg, W2, v0, invTau, invSigma2 ); 
	
	if( commonpi ){
		if( all(beta1pi == 0.0) ){
			priorInclusion.fill( arma::accu(incProbs) / double(JD) ); 
		}else if( all(beta1pi > 0.0) ){
			double sumtmp; 
			betastar1.fill( beta1pi(0) + arma::accu( incProbs ) ); 
			betastar2.fill( beta2pi(0) + arma::accu( 1.0 - incProbs) ); 
			sumtmp = betastar1(0) + betastar2(0); 
			priorInclusion.fill( betastar1(0) / (sumtmp) );
		}else{} 			// Don't update probabilities if priorInclusion is fixed    
	}else{
		int d; 
		if( all(beta1pi == 0.0) ){
			for( d=0; d<D; d++ ){
				priorInclusion(d) = ( arma::accu(incProbs.col(d)) / double(J) ); 	
			}
		}else if( all(beta1pi > 0.0) ){
			double sumtmp; 
			for( d=0; d<D; d++ ){
				betastar1(d) = beta1pi(d) + arma::accu( incProbs.col(d) ); 
				betastar2(d) = beta2pi(d) + arma::accu( 1.0 - incProbs.col(d));
				sumtmp = betastar1(d) + betastar2(d); 		
				priorInclusion(d) = ( betastar1(d) / (sumtmp) );			
			}			
		}else{}				// Don't update probabilities if priorInclusion is fixed  	
	}
}







