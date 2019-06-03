#ifndef VBALGORITHM_H
#define VBALGORITHM_H

void vbalgorithm( const arma::mat& X, int D, int I, int J, int maxIter, double tolerance, 
                  bool updatetau, std::string priorvar,
                  arma::vec alphatau, arma::vec betatau, arma::vec gammatau, arma::vec deltatau, bool SVS, arma::vec priorInclusion, arma::vec beta1pi, arma::vec beta2pi, double v0, 
				  std::string hypertype, bool commonpi, int JD, 
				  arma::mat invTau, double qz, bool scaleprior, Rcpp::List hpdis, int &it, double sigma2, double aPostSigma, double denomX, double normX, double elbo, 
				  double hW, double &finalElbo, double EWtauW, bool &converged, 
				  arma::vec betastar1, arma::vec betastar2, arma::mat muW, arma::mat muP, 
				  arma::mat W2, const arma::mat& XTX, arma::mat incProbs, arma::mat deltataupost,
                  double &globalElbo, arma::mat &globalMuW, arma::mat &globalMuP,
                  arma::mat &globaltau, double &globalSigma2, 
                  arma::mat &globalbetatau, bool &globalconverged, arma::vec &elbovals, 
                  Rcpp::List &globalhpdi, arma::vec &globalPriorInc, 
                  arma::mat &globalIncPr, bool globalvar, bool hpdi, bool verbose ); 

#endif