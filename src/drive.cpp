#include "drive.h"

#ifndef CCODE
using Rcpp::Rcerr;
using Rcpp::Rcout;
#else
#define Rcout std::cout
#define Rcerr std::cerr
#endif

#ifdef _OPENMP
extern omp_lock_t RNGlock; /*defined in global.h*/
#endif
/*extern std::vector<std::mt19937_64> rng;*/

// not yet set up openmp
#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]

using Utils::Chain_Data;

// main function
// (i) import data and parameters; (ii) MCMC algorithm; (iii) export estimates
Rcpp::List drive(const std::string &dataFile, const unsigned int p, const unsigned int q, const std::string &hyperParFile, const std::string &outFilePath,
                 const arma::vec ini_beta, const arma::vec ini_tauSq, const arma::vec ini_h, const arma::uvec groupInd,
                 const unsigned int nIter, const unsigned int nChains, const unsigned int thin, bool rw)
{

// set random seed
#ifdef _OPENMP
    // NOT yet/easily implement for updateRP().
    // Through a simple test, '#pragma omp parallel for' slowed down the computation of other small functions
    /*omp_set_nested( 0 );
    omp_set_num_threads( 1 );*/
    // this slows down running on Linux machines
    omp_init_lock(&RNGlock); // init RNG lock for the parallel part
    // arma::arma_rng::set_seed(123);
#endif

    // ###########################################################
    // ## Read Arguments and Data
    // ###########################################################

    // Rcout << std::fixed << std::setprecision(11);

    Chain_Data chainData; // this initialises the pointers and the strings to ""
    chainData.outFilePath = outFilePath;

    PSBC myPSBC;

    // read Data and format into usables
    // Rcout << "Reading input files ... " <<  "\n";

    Utils::formatData(dataFile, chainData.survData);
    Utils::readHyperPar(hyperParFile, chainData);

    // declare all the data-related variables
    double eta0 = chainData.eta0;
    double kappa0 = chainData.kappa0;
    double c0 = chainData.c0;
    double r = chainData.r;
    double delta = chainData.delta;

    arma::vec dataTime = chainData.survData.data->col(0);
    arma::vec dataDi = chainData.survData.data->col(1);
    arma::vec s0 = {0.}; // event times that are actually observed
    s0(0) = 2. * max(dataTime) - max(dataTime.elem(arma::find(dataTime != max(dataTime))));
    arma::vec s = arma::join_cols(arma::sort(arma::unique(dataTime.elem(arma::find(dataDi == 1.)))), s0);
    unsigned int J = s.n_elem;
    arma::uvec groupNo = arma::unique(groupInd);
    unsigned int K = groupNo.n_elem;

    ind_r_d = ind_r = ind_d = arma::zeros<arma::mat>(chainData.survData.data->n_rows, J);
    d = arma::sum(ind_d.t(), 1);
    myPSBC.settingInterval_cpp(dataTime, dataDi, s, J, ind_d, ind_r, ind_r_d, d);

    arma::vec beta_prop_me;
    be = beta_prop_me = ini_beta;
    arma::mat beta_p = arma::zeros<arma::mat>((int)(nIter / thin) + 1, p + q);
    beta_p.row(0) = ini_beta.t();

    double lambdaSq;
    lambdaSq = chainData.lambdaSq;
    double sigmaSq;
    sigmaSq = chainData.sigmaSq;
    arma::vec tauSq = ini_tauSq;
    h = ini_h;

    // tausq: only for genomic variables
    arma::vec tauSq_exp = arma::zeros<arma::vec>(p);
    for (unsigned int i = 0; i < K; ++i)
    {
        // groupNoLocation = arma::find( groupInd == groupNo(i) );
        tauSq_exp.elem(arma::find(groupInd == groupNo(i))).fill(tauSq(i));
    }

    double beta_prop_sd = sqrt(chainData.beta_prop_var);
    double beta_clin_sd = sqrt(chainData.beta_clin_var);
    arma::vec sd_bePart2 = arma::zeros<arma::vec>(q);
    sd_bePart2.fill(beta_clin_sd);
    arma::vec sd_be = arma::join_cols(sqrt(sigmaSq * tauSq_exp), sd_bePart2);

    // arma::uvec xId = arma::linspace<arma::uvec>(2, 1+p+q, p+q)
    arma::mat x = chainData.survData.data->cols(arma::linspace<arma::uvec>(2, 1 + p + q, p + q));
    arma::vec xbeta = x * ini_beta;

    arma::vec be_normSq = arma::zeros<arma::vec>(K);
    for (unsigned int i = 0; i < K; ++i)
        be_normSq(i) = arma::accu(ini_beta.elem(arma::find(groupInd == groupNo(i))) % ini_beta.elem(arma::find(groupInd == groupNo(i))));

    if (!any(be_normSq))
        be_normSq(0) = 0.1;

    arma::vec alpha0 = arma::zeros<arma::vec>(1 + J);
    // H_star(j) = eta0 * arma::pow( s, kappa0 );
    alpha0.subvec(1, J) = c0 * eta0 * arma::pow(s, kappa0);
    arma::vec hPriorSh = arma::diff(alpha0);

    // save mcmc intermediate estimates
    arma::mat h_p = arma::zeros<arma::mat>((int)(nIter / thin) + 1, h.n_elem);
    h_p.row(0) = h.t();
    arma::mat tauSq_p = arma::zeros<arma::mat>((int)(nIter / thin) + 1, tauSq.n_elem);
    tauSq_p.row(0) = tauSq.t();

    arma::vec sigmaSq_p = arma::zeros<arma::vec>((int)(nIter / thin) + 1);
    arma::vec lambdaSq_p = sigmaSq_p;
    sampleRPc_accept = arma::zeros<arma::uvec>(q);
    sampleRPg_accept = arma::zeros<arma::uvec>(p);

    // ###########################################################
    // ## MCMC loop
    // ###########################################################

    // unsigned int tick = 100; // how many iter for each print?
    unsigned int j = 0; // count thinned results
    const int cTotalLength = 50;
    Rcout << "Running MCMC iterations ...\n";

    for (unsigned int M = 0; M < nIter; ++M)
    {

        if (M % 10 == 0 || M == (nIter - 1))
            Rcout << "\r[" <<                                               //'\r' aka carriage return should move printer's cursor back at the beginning of the current line
                std::string(cTotalLength * (M + 1.) / nIter, '#') <<        // printing filled part
                std::string(cTotalLength * (1. - (M + 1.) / nIter), '-') << // printing empty part
                "] " << int((M + 1.) / nIter * 100.0) << "%\r";             // printing percentage

        // Updating regression coefficients and hyperparameters

        if (q > 0)
            myPSBC.updateRP_clinical_cpp(p, q, x, ind_r, ind_d, ind_r_d, J, beta_prop_me, beta_prop_sd, xbeta, be, h, sd_be, sampleRPc_accept);

        if (rw)
        {
            myPSBC.updateRP_genomic_rw_cpp(p, x, ind_r, ind_d, ind_r_d, J, beta_prop_me, beta_prop_sd, xbeta, be, h, sd_be, sampleRPg_accept);
        }
        else
        {
            myPSBC.updateRP_genomic_cpp(p, x, ind_r, ind_d, ind_r_d, J, xbeta, be, h, sd_be, sampleRPg_accept);
        }

        //        if( q > 0 )
        //        {
        for (unsigned int i = 0; i < K; ++i)
            be_normSq(i) = arma::accu(be.elem(arma::find(groupInd == groupNo(i))) % be.elem(arma::find(groupInd == groupNo(i))));
        //        } else {
        //            be_normSq = be % be;
        //        }

        h = myPSBC.updateBH_cpp(ind_r_d, hPriorSh, d, c0, J, xbeta);

        tauSq = myPSBC.updateTau_GL_cpp(lambdaSq, sigmaSq, be_normSq);

        sigmaSq = myPSBC.updateSigma_GL_cpp(p, be_normSq, tauSq);
        lambdaSq = myPSBC.updateLambda_GL_cpp(p, K, r, delta, tauSq);

        //        if( q > 0 )
        //        {
        for (unsigned int i = 0; i < K; ++i)
            tauSq_exp.elem(arma::find(groupInd == groupNo(i))).fill(tauSq(i));

        sd_bePart2.fill(beta_clin_sd);
        sd_be = arma::join_cols(sqrt(sigmaSq * tauSq_exp), sd_bePart2);
        //        } else {
        //            tauSq_exp = tauSq;
        //            sd_be = sqrt(sigmaSq * tauSq_exp);
        //        }

        // Save all results
        if (M % thin == 0)
        {
            ++j;
            beta_p.row(j) = be.t();
            h_p.row(j) = h.t();
            tauSq_p.row(j) = tauSq.t();
            sigmaSq_p(j) = sigmaSq;
            lambdaSq_p(j) = lambdaSq;
        }

        // adaptve jumping rule
        if (j > 20)
            for (unsigned int jj = 0; jj < p + q; ++jj)
                if (ini_beta(jj) == beta_p(j - 20, jj))
                    beta_prop_me(jj) = be(jj);
    }

    arma::uvec accept_rate = join_cols(sampleRPg_accept, sampleRPc_accept); // / (double)(nIter);

    // Exit
    Rcout << "\nDONE, exiting! \n";

    return Rcpp::List::create(Rcpp::Named("beta.p") = beta_p,
                              Rcpp::Named("h.p") = h_p,
                              Rcpp::Named("tauSq.p") = tauSq_p,
                              Rcpp::Named("sigmaSq.p") = sigmaSq_p,
                              Rcpp::Named("lambdaSq.p") = lambdaSq_p,
                              Rcpp::Named("accept.rate") = accept_rate);
}

// [[Rcpp::export]]
Rcpp::List psbcSpeedUp_internal(const std::string &dataFile, const unsigned int p, const unsigned int q, const std::string &hyperParFile, const std::string &outFilePath,
                                const arma::vec ini_beta, const arma::vec ini_tauSq, const arma::vec ini_h, const arma::uvec groupInd, const unsigned int nIter, const unsigned int nChains, const unsigned int thin, bool rw)
{
    // int status {1};
    Rcpp::List beta_mcmc;
    try
    {
        // status = drive(...);
        beta_mcmc = drive(dataFile, p, q, hyperParFile, outFilePath, ini_beta, ini_tauSq, ini_h, groupInd, nIter, nChains, thin, rw);
    }
    catch (const std::exception &e)
    {
        Rcout << e.what() << '\n'; // we can use Rcerr here because we're reaching here from R for sure
    }

    // return status;
    return beta_mcmc;
}
