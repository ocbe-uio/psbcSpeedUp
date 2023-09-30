#ifndef DRIVEPSBC
#define DRIVEPSBC

#include "global.h"
#include "utils.h"

#include <vector>
#include <sstream>
#include <random>

arma::mat vecToMat( const arma::vec, const int );
arma::mat matProdVec( const arma::vec, const arma::vec& );
void settingInterval_cpp( const arma::vec, const arma::vec, const arma::vec, unsigned int, arma::mat&, arma::mat&, arma::mat&, arma::vec& );
double updateLambda_GL_cpp( int, unsigned int, double, double, arma::vec );
double updateSigma_GL_cpp( int, arma::vec, arma::vec );
arma::vec rinvgauss( arma::vec, double );
arma::vec updateTau_GL_cpp( double, double, arma::vec );
arma::vec updateBH_cpp( arma::mat&, arma::mat, arma::mat&, int, unsigned int, arma::vec );
void updateRP_clinical_cpp( int, int, const arma::mat, arma::mat&, arma::mat&, arma::mat&, unsigned int, arma::vec, double, arma::vec&, arma::vec&, arma::vec&, arma::vec, arma::uvec& );
void updateRP_genomic_cpp( int, const arma::mat, arma::mat&, arma::mat&, arma::mat&, unsigned int, arma::vec&, arma::vec&, arma::vec&, arma::vec, arma::uvec& );
void updateRP_genomic_rw_cpp( int, const arma::mat, arma::mat&, arma::mat&, arma::mat&, unsigned int, arma::vec, double, arma::vec&, arma::vec&, arma::vec&, arma::vec, arma::uvec& );

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

Rcpp::List drive( const std::string& dataFile, const int p, const int q, const std::string& hyperParFile, const std::string& outFilePath,
                const arma::vec ini_beta, const arma::vec ini_tauSq, const arma::vec ini_h, const arma::uvec groupInd,
          unsigned int nIter, unsigned int nChains, unsigned int thin, bool rw);

Rcpp::List psbcSpeedUp_internal( const std::string& dataFile, const int p, const int q, const std::string& hyperParFile, //const std::string& initialParFile,
                         const std::string& outFilePath,
                                const arma::vec ini_beta, const arma::vec ini_tauSq, const arma::vec ini_h, const arma::uvec groupInd,
                         unsigned int nIter=10, unsigned int nChains=1, unsigned int thin=1, bool rw=false );
#endif
