#ifndef DRIVEPSBC
#define DRIVEPSBC

#include "global.h"
#include "psbc.h"

Rcpp::List drive(
    const arma::vec& datTime, 
    const arma::uvec& datEvent, 
    const arma::mat& x,
    const unsigned int p, 
    const unsigned int q, 
    Rcpp::List &hyperParFile,
    arma::vec beta, 
    arma::vec tauSq, 
    arma::vec h, 
    const arma::uvec groupInd, 
    const unsigned int nIter, 
    const unsigned int nChains, 
    const unsigned int thin, 
    bool rw
);

#endif
