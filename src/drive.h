#ifndef DRIVEPSBC
#define DRIVEPSBC

#include "global.h"
#include "psbc.h"

Rcpp::List drive( 
    const arma::vec& datTime, const arma::uvec& datEvent, const arma::mat& x,
    const unsigned int p, const unsigned int q, Rcpp::List &hyperParFile, 
                 arma::vec ini_beta, arma::vec ini_tauSq, arma::vec ini_h, const arma::uvec groupInd,
                 const unsigned int nIter, const unsigned int nChains, const unsigned int thin, bool rw);

Rcpp::List psbcSpeedUp_internal( 
    const arma::vec& datTime, const arma::uvec& datEvent, const arma::mat& x,
    const unsigned int p, const unsigned int q, Rcpp::List &hyperParFile, 
                                arma::vec ini_beta, arma::vec ini_tauSq, arma::vec ini_h, const arma::uvec groupInd,
                                const unsigned int nIter = 10, const unsigned int nChains = 1, const unsigned int thin = 1, bool rw = false);


#endif
