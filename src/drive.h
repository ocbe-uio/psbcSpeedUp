#ifndef DRIVEPSBC
#define DRIVEPSBC

#include "global.h"
#include "utils.h"
#include "psbc.h"

Rcpp::List drive(const std::string &dataFile, const unsigned int p, const unsigned int q, Rcpp::List &hyperParFile, 
                 const arma::vec ini_beta, const arma::vec ini_tauSq, const arma::vec ini_h, const arma::uvec groupInd,
                 const unsigned int nIter, const unsigned int nChains, const unsigned int thin, bool rw);

Rcpp::List psbcSpeedUp_internal(const std::string &dataFile, const unsigned int p, const unsigned int q, Rcpp::List &hyperParFile, 
                                const arma::vec ini_beta, const arma::vec ini_tauSq, const arma::vec ini_h, const arma::uvec groupInd,
                                const unsigned int nIter = 10, const unsigned int nChains = 1, const unsigned int thin = 1, bool rw = false);


#endif
