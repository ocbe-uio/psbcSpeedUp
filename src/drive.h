#ifndef DRIVEPSBC
#define DRIVEPSBC

#include "global.h"
#include "utils.h"
#include "psbc.h"

Rcpp::List drive( const std::string& dataFile, const unsigned int p, const unsigned int q, const std::string& hyperParFile, const std::string& outFilePath,
                const arma::vec ini_beta, const arma::vec ini_tauSq, const arma::vec ini_h, const arma::uvec groupInd,
                 const unsigned int nIter, const unsigned int nChains, const unsigned int thin, bool rw);

Rcpp::List psbcSpeedUp_internal( const std::string& dataFile, const unsigned int p, const unsigned int q, const std::string& hyperParFile, //const std::string& initialParFile,
                         const std::string& outFilePath,
                                const arma::vec ini_beta, const arma::vec ini_tauSq, const arma::vec ini_h, const arma::uvec groupInd,
                                const unsigned int nIter=10, const unsigned int nChains=1, const unsigned int thin=1, bool rw=false );


// define global variables
arma::mat ind_r;
arma::mat ind_d;
arma::mat ind_r_d;

arma::vec xbeta;
arma::vec be;
arma::uvec sampleRPc_accept, sampleRPg_accept;
arma::vec d; // used to update "h"
arma::vec h;



//inline arma::mat& h_exp_xbeta_mat = {0.};
//inline arma::mat& h_exp_xbeta_prop_mat = {0.};

#endif
