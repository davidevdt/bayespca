#ifndef UPDATEELBO_H
#define UPDATEELBO_H

void updateElbo( double &elbo, double sigma2, bool globalvar, std::string priorvar,
                 arma::mat W2, arma::vec betatau, arma::vec alphatau,
                 int J, int D,
                 arma::mat Tau, int I, double denomX, double hW, int JD );
#endif
