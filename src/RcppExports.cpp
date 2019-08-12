#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mainBayesPCA
Rcpp::List mainBayesPCA(const arma::mat& X, int D, int I, int J, int nstart, int maxIter, double tolerance, bool svdStart, bool verbose, bool updatetau, std::string priorvar, arma::vec alphatau, arma::vec betatau, arma::vec gammatau, arma::vec deltatau, bool SVS, arma::vec priorInclusion, arma::vec beta1pi, arma::vec beta2pi, double v0, bool commonpi, int JD, arma::mat invTau, double qz, bool scaleprior, std::string hypertype, bool globalvar, bool hpdi);
RcppExport SEXP _bayespca_mainBayesPCA(SEXP XSEXP, SEXP DSEXP, SEXP ISEXP, SEXP JSEXP, SEXP nstartSEXP, SEXP maxIterSEXP, SEXP toleranceSEXP, SEXP svdStartSEXP, SEXP verboseSEXP, SEXP updatetauSEXP, SEXP priorvarSEXP, SEXP alphatauSEXP, SEXP betatauSEXP, SEXP gammatauSEXP, SEXP deltatauSEXP, SEXP SVSSEXP, SEXP priorInclusionSEXP, SEXP beta1piSEXP, SEXP beta2piSEXP, SEXP v0SEXP, SEXP commonpiSEXP, SEXP JDSEXP, SEXP invTauSEXP, SEXP qzSEXP, SEXP scalepriorSEXP, SEXP hypertypeSEXP, SEXP globalvarSEXP, SEXP hpdiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type I(ISEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type nstart(nstartSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< bool >::type svdStart(svdStartSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type updatetau(updatetauSEXP);
    Rcpp::traits::input_parameter< std::string >::type priorvar(priorvarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alphatau(alphatauSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type betatau(betatauSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gammatau(gammatauSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type deltatau(deltatauSEXP);
    Rcpp::traits::input_parameter< bool >::type SVS(SVSSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type priorInclusion(priorInclusionSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta1pi(beta1piSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta2pi(beta2piSEXP);
    Rcpp::traits::input_parameter< double >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< bool >::type commonpi(commonpiSEXP);
    Rcpp::traits::input_parameter< int >::type JD(JDSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invTau(invTauSEXP);
    Rcpp::traits::input_parameter< double >::type qz(qzSEXP);
    Rcpp::traits::input_parameter< bool >::type scaleprior(scalepriorSEXP);
    Rcpp::traits::input_parameter< std::string >::type hypertype(hypertypeSEXP);
    Rcpp::traits::input_parameter< bool >::type globalvar(globalvarSEXP);
    Rcpp::traits::input_parameter< bool >::type hpdi(hpdiSEXP);
    rcpp_result_gen = Rcpp::wrap(mainBayesPCA(X, D, I, J, nstart, maxIter, tolerance, svdStart, verbose, updatetau, priorvar, alphatau, betatau, gammatau, deltatau, SVS, priorInclusion, beta1pi, beta2pi, v0, commonpi, JD, invTau, qz, scaleprior, hypertype, globalvar, hpdi));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayespca_mainBayesPCA", (DL_FUNC) &_bayespca_mainBayesPCA, 28},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayespca(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
