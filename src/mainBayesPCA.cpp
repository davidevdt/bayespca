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
                        arma::vec gammatau, arma::vec deltatau, bool SVS,
                        arma::vec priorInclusion, arma::vec beta1pi, arma::vec beta2pi,
                        double v0, bool commonpi,
                        int JD, arma::mat invTau, double qz, bool scaleprior, 
                        std::string hypertype, bool globalvar, bool hpdi 
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
  
	aPostSigma = ( double(J) * double(I)) - double(J); 
  
	if( scaleprior ){
		aPostSigma += double(JD);
	}
  
	double denomX = 0.0;  
	double normX; 
	double elbo = 0.0; 
	double hW = 0.0; 
	double finalElbo = 0.0; 
	double EWtauW = 0.0; 
  
	// Initializations: Booleans
	bool converged = FALSE; 
  
	// Initializations: Vectors
	arma::vec betastar1; 
	arma::vec betastar2;  

	// Initializations : Matrices
	arma::mat muW(J, D);
	arma::mat muP(J, D); 
	arma::mat W2(J, D);
	arma::mat XTX(J, J); 
	arma::mat incProbs; 
  
  
	arma::mat deltataupost;

	if( globalvar == FALSE ){
		if( (priorvar != "invgamma" ) || (any( gammatau < 0.0 ))  ){
			deltataupost.set_size(1,1); 
		}else{
			if( hypertype == "common" ){
				deltataupost.set_size(1,1);
			}else if( hypertype == "component" ){
				deltataupost.set_size(1,D); 
			}else{
				deltataupost.set_size(J,D); 
			}
		}
	}else{
		if( priorvar == "invgamma" && all( gammatau > 0.0 ) ){
			deltataupost.set_size(1,D);
			deltataupost.row(0) = deltatau.t(); 
		}else{  
			deltataupost.set_size(1,1);
		}
	}


	if( SVS ){
		incProbs.set_size( J,D ); 

		if( all(beta1pi > 0.0) ){
			betastar1.set_size( D ); 
			betastar2.set_size( D );
		}else{
			betastar1.set_size( 1 ); 
			betastar2.set_size( 1 );			
		}
			
		if( commonpi ){
			incProbs.fill( priorInclusion(0) );

			betastar1.fill( beta1pi(0) + arma::accu( incProbs ) ); 
			betastar2.fill( beta2pi(0) + arma::accu( 1.0 - incProbs ) ); 
			
			
		}else{

			incProbs.each_row() = priorInclusion.t(); 
			
			if( all(beta1pi > 0.0) ){

				betastar1 = beta1pi + arma::sum( incProbs, 0 ).t(); 
				betastar2 = beta2pi + arma::sum( 1.0 - incProbs, 0 ).t();
			}
		}	
	}else{
		incProbs.set_size( 1,1 ); 
		betastar1.set_size(1); 
		betastar2.set_size(1);  
	}

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
	arma::mat globalbetatau;
	bool globalconverged = FALSE;
	arma::vec elbovals; 
	Rcpp::List globalhpdi; 
	arma::vec globalPriorInc(D); 
	arma::mat globalIncPr; 

	
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
					priorvar, alphatau, betatau, gammatau, deltatau, 
					SVS, priorInclusion, beta1pi, beta2pi, v0, 
					hypertype, commonpi, JD, invTau, qz, scaleprior, hpdis, 
					it, sigma2, aPostSigma, denomX, normX, elbo, hW,
					finalElbo, EWtauW, converged, betastar1, betastar2, muW, muP, W2, XTX, 	incProbs, deltataupost, globalElbo, globalMuW, globalMuP, globaltau, 
					globalSigma2, globalbetatau, globalconverged, elbovals, 
					globalhpdi, globalPriorInc, globalIncPr, 
					globalvar, hpdi, verbose );			 
    
		
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
	ret["globalbetatau"] = globalbetatau; 
	ret["globalPriorInc"] = globalPriorInc; 
	ret["globalIncPr"] = globalIncPr;
 
  
	return ret; 
  
} 
