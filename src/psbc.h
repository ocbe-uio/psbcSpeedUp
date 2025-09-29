#ifndef PSBC_H
#define PSBC_H

#include <vector>
#include <sstream>
#include <random>

#ifdef CCODE
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif

class PSBC
{

public:
    PSBC();
    ~PSBC();

    static void settingInterval_cpp(const arma::vec& y, const arma::uvec& delta_, const arma::vec s_, const unsigned int J_, arma::mat &ind_d_, arma::mat &ind_r_, arma::mat &ind_r_d_, arma::vec &d_);
    static double updateLambda_GL_cpp(const unsigned int, const unsigned int, double, double, arma::vec);
    static double updateSigma_GL_cpp(const unsigned int, arma::vec, arma::vec);
    static arma::vec updateTau_GL_cpp(double, double, arma::vec);
    static arma::vec updateBH_cpp(arma::mat &, arma::vec, arma::vec &, double, const unsigned int, arma::vec);
    
    static void updateRP_clinical_cpp(const unsigned int p, const unsigned int q, const arma::mat& x_, arma::mat &ind_r_, arma::mat &ind_d_, arma::mat &ind_r_d_, const unsigned int J_, arma::vec beta_prop_me_, double beta_prop_sd, arma::vec &xbeta_, arma::vec &be_, arma::vec &h_, arma::vec sd_be_, arma::uvec &sampleRPc_accept_);
    static void updateRP_genomic_cpp(const unsigned int p, const arma::mat& x_, arma::mat &ind_r_, arma::mat &ind_d_, arma::mat &ind_r_d_, const unsigned int J_, arma::vec &xbeta_, arma::vec &be_, arma::vec &h_, arma::vec sd_be_, arma::uvec &sampleRPg_accept_);
    static void updateRP_genomic_rw_cpp(const unsigned int p, const arma::mat& x_, arma::mat &ind_r_, arma::mat &ind_d_, arma::mat &ind_r_d_, const unsigned int J_, arma::vec beta_prop_me_, double beta_prop_sd, arma::vec &xbeta_, arma::vec &be_, arma::vec &h_, arma::vec sd_be_, arma::uvec &sampleRPg_accept_);

private:

    static arma::mat matProdVec(const arma::mat, const arma::vec);
    static arma::vec sumMatProdVec(const arma::mat, const arma::vec);
    static arma::vec rinvgauss(arma::vec, double);
    
};

#endif
