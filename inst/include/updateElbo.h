#ifndef UPDATEELBO_H
#define UPDATEELBO_H

void updateElbo( double &elbo, double sigma2, bool globalvar, std::string priorvar,
<<<<<<< HEAD
                 arma::mat W2, bool SVS, arma::mat incProbs, double v0, arma::vec gammatau,
                 arma::vec betatau, arma::mat deltataupost, arma::vec alphatau,
                 arma::vec deltatau, std::string hypertype, int J, int D,
                 arma::mat Tau, int I, double denomX, double hW, int JD, 
                 bool commonpi, arma::vec beta1pi, arma::vec beta2pi, arma::vec priorInclusion,
                 arma::vec betastar1, arma::vec betastar2, bool scaleprior );
=======
                 arma::mat W2, arma::vec betatau, arma::vec alphatau,
                 int J, int D,
                 arma::mat Tau, int I, double denomX, double hW, int JD );
>>>>>>> 50009e97c685ef8e94bbfdb6fc3a466f64df3285
#endif
