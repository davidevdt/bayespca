#ifndef UPDATEEXPECTED_H
#define UPDATEEXPECTED_H

void updateExpected( arma::mat &muW, arma::mat &W2, double &sigma2, arma::mat &muP, int D,
										 int J, int I, double &hW, const arma::mat& X, const arma::mat& XTX,
										 arma::mat &Tau, bool hpdi, Rcpp::List &hpdis, double qz, double &denomX,
										 double aPostSigma, double normX, std::string priorvar,
										 bool updatetau, bool globalvar, int JD, 
										 arma::vec alphatau, arma::vec betatau );

#endif
