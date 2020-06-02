#ifndef VBALGORITHM_H
#define VBALGORITHM_H

void vbalgorithm( const arma::mat& X, int D, int I, int J, int maxIter, double tolerance,
                  bool updatetau, std::string priorvar, arma::vec alphatau, arma::vec betatau,
                  int JD, arma::mat Tau, double qz, Rcpp::List hpdis, int &it, double sigma2,
                  double aPostSigma, double denomX, double normX, double elbo,
				          double hW, double &finalElbo, bool &converged,
                  arma::mat muW, arma::mat muP, arma::mat W2, const arma::mat& XTX,
                  double &globalElbo, arma::mat &globalMuW, arma::mat &globalMuP,
                  arma::mat &globaltau, double &globalSigma2,
                  bool &globalconverged, arma::vec &elbovals,
                  Rcpp::List &globalhpdi, bool globalvar, bool hpdi, bool verbose );

#endif
