#ifndef UPDATEEXPECTED_H
#define UPDATEEXPECTED_H

void updateExpected( arma::mat &muW, arma::mat &W2, double &sigma2, arma::mat &muP, int D,
<<<<<<< HEAD
					 int J, int I, double &hW, const arma::mat& X, const arma::mat& XTX,
					 arma::mat &Tau, bool hpdi,
					 Rcpp::List &hpdis, double qz, bool scaleprior, double &EWtauW, double v0,
					 arma::mat &incProbs, double SVS, double &denomX,
					 double aPostSigma, double normX,
					 std::string priorvar,
					 bool updatetau, bool globalvar, arma::vec gammatau, arma::vec alphatau,
					 arma::vec &betatau, std::string hypertype, arma::mat &deltataupost,
					 arma::vec deltatau, int JD, arma::vec &priorInclusion,
					 arma::vec &betastar1, arma::vec &betastar2, bool commonpi,
					 arma::vec beta1pi, arma::vec beta2pi
					);
=======
										 int J, int I, double &hW, const arma::mat& X, const arma::mat& XTX,
										 arma::mat &Tau, bool hpdi, Rcpp::List &hpdis, double qz, double &denomX,
										 double aPostSigma, double normX, std::string priorvar,
										 bool updatetau, bool globalvar, int JD, 
										 arma::vec alphatau, arma::vec betatau );
>>>>>>> 50009e97c685ef8e94bbfdb6fc3a466f64df3285

#endif
