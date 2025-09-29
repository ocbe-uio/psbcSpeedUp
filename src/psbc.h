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

    static void settingInterval_cpp(const arma::vec, const arma::vec, const arma::vec, const unsigned int, arma::mat &, arma::mat &, arma::mat &, arma::vec &);
    static double updateLambda_GL_cpp(const unsigned int, const unsigned int, double, double, arma::vec);
    static double updateSigma_GL_cpp(const unsigned int, arma::vec, arma::vec);
    static arma::vec updateTau_GL_cpp(double, double, arma::vec);
    static arma::vec updateBH_cpp(arma::mat &, arma::vec, arma::vec &, double, const unsigned int, arma::vec);
    static void updateRP_clinical_cpp(const unsigned int, const unsigned int, const arma::mat, arma::mat &, arma::mat &, arma::mat &, const unsigned int, arma::vec, double, arma::vec &, arma::vec &, arma::vec &, arma::vec, arma::uvec &);
    static void updateRP_genomic_cpp(const unsigned int, const arma::mat, arma::mat &, arma::mat &, arma::mat &, const unsigned int, arma::vec &, arma::vec &, arma::vec &, arma::vec, arma::uvec &);
    static void updateRP_genomic_rw_cpp(const unsigned int, const arma::mat, arma::mat &, arma::mat &, arma::mat &, const unsigned int, arma::vec, double, arma::vec &, arma::vec &, arma::vec &, arma::vec, arma::uvec &);

private:

    static arma::mat matProdVec(const arma::mat, const arma::vec);
    static arma::vec sumMatProdVec(const arma::mat, const arma::vec);
    static arma::vec rinvgauss(arma::vec, double);
    
};

#endif
