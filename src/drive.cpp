#include "drive.h"

#ifdef _OPENMP
extern omp_lock_t RNGlock; /*defined in global.h*/
#endif
// extern std::vector<std::mt19937_64> rng;

// TODO: not yet set up openmp
#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]


// main function
// (i) import data and parameters; (ii) MCMC algorithm; (iii) export estimates

// [[Rcpp::export]]
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
    bool rw)
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

    // get hyperparameters
    double hyperPar_eta0 = Rcpp::as<double>(hyperParFile["eta0"]);
    double hyperPar_kappa0 = Rcpp::as<double>(hyperParFile["kappa0"]);
    double hyperPar_c0 = Rcpp::as<double>(hyperParFile["c0"]);
    double hyperPar_r = Rcpp::as<double>(hyperParFile["r"]);
    double hyperPar_delta = Rcpp::as<double>(hyperParFile["delta"]);
    double hyperPar_beta_prop_sd = sqrt(Rcpp::as<double>(hyperParFile["beta_prop_var"]));
    double hyperPar_beta_clin_sd = sqrt(Rcpp::as<double>(hyperParFile["beta_clin_var"]));
    double hyperPar_lambdaSq = Rcpp::as<double>(hyperParFile["lambdaSq"]);
    // double rate = Rcpp::as<double>(hyperPar["rate"]);
    double hyperPar_sigmaSq = Rcpp::as<double>(hyperParFile["sigmaSq"]);
    hyperParFile = Rcpp::List();  // Clear it by creating a new empty List

    // arma::vec datTime = chainData.survData.data->col(0);
    // arma::uvec datEvent = arma::conv_to<arma::ivec>(chainData.survData.data->col(1));
    arma::vec s0 = {0.}; // event times that are actually observed
    s0(0) = 2. * max(datTime) - max(datTime.elem(arma::find(datTime != max(datTime))));
    arma::vec s = arma::join_cols(arma::sort(arma::unique(datTime.elem(arma::find(datEvent == 1)))), s0);
    unsigned int J = s.n_elem;
    arma::uvec groupNo = arma::unique(groupInd);
    unsigned int K = groupNo.n_elem;

    arma::mat ind_r_d = arma::zeros<arma::mat>(datEvent.n_elem, J);
    arma::mat ind_r = arma::zeros<arma::mat>(datEvent.n_elem, J);
    arma::mat ind_d = arma::zeros<arma::mat>(datEvent.n_elem, J);
    arma::vec d = arma::sum(ind_d.t(), 1);
    PSBC::settingInterval_cpp(datTime, datEvent, s, J, ind_d, ind_r, ind_r_d, d);

    arma::vec beta_prop_me = beta;
    arma::mat beta_p = arma::zeros<arma::mat>((int)(nIter / thin) + 1, p + q);
    beta_p.row(0) = beta.t();

    // tausq: only for genomic variables
    arma::vec tauSq_exp = arma::zeros<arma::vec>(p);
    for (unsigned int i = 0; i < K; ++i)
    {
        // groupNoLocation = arma::find( groupInd == groupNo(i) );
        tauSq_exp.elem(arma::find(groupInd == groupNo(i))).fill(tauSq(i));
    }

    arma::vec sd_bePart2 = arma::zeros<arma::vec>(q);
    sd_bePart2.fill(hyperPar_beta_clin_sd);
    arma::vec sd_be = arma::join_cols(sqrt(hyperPar_sigmaSq * tauSq_exp), sd_bePart2);

    // arma::uvec xId = arma::linspace<arma::uvec>(2, 1+p+q, p+q)
    // arma::mat x = chainData.survData.data->cols(arma::linspace<arma::uvec>(2, 1 + p + q, p + q));
    arma::vec xbeta = x * beta;

    arma::vec be_normSq = arma::zeros<arma::vec>(K);
    for (unsigned int i = 0; i < K; ++i)
        be_normSq(i) = arma::accu(beta.elem(arma::find(groupInd == groupNo(i))) % beta.elem(arma::find(groupInd == groupNo(i))));

    if (!any(be_normSq))
        be_normSq(0) = 0.1;

    arma::vec alpha0 = arma::zeros<arma::vec>(1 + J);
    // H_star(j) = eta0 * arma::pow( s, kappa0 );
    alpha0.subvec(1, J) = hyperPar_c0 * hyperPar_eta0 * arma::pow(s, hyperPar_kappa0);
    arma::vec hPriorSh = arma::diff(alpha0);

    // save mcmc intermediate estimates
    arma::mat h_p = arma::zeros<arma::mat>((int)(nIter / thin) + 1, h.n_elem);
    h_p.row(0) = h.t();
    arma::mat tauSq_p = arma::zeros<arma::mat>((int)(nIter / thin) + 1, tauSq.n_elem);
    tauSq_p.row(0) = tauSq.t();

    arma::vec sigmaSq_p = arma::zeros<arma::vec>((int)(nIter / thin) + 1);
    arma::vec lambdaSq_p = sigmaSq_p;
    arma::uvec sampleRPc_accept = arma::zeros<arma::uvec>(q);
    arma::uvec sampleRPg_accept = arma::zeros<arma::uvec>(p);

    // ###########################################################
    // ## MCMC loop
    // ###########################################################

    // unsigned int tick = 100; // how many iter for each print?
    unsigned int j = 0; // count thinned results
    const int cTotalLength = 50;
    Rcpp::Rcout << "Running MCMC iterations ...\n";

    for (unsigned int M = 0; M < nIter; ++M)
    {

        if (M % 10 == 0 || M == (nIter - 1))
            Rcpp::Rcout << "\r[" <<                                               //'\r' aka carriage return should move printer's cursor back at the beginning of the current line
                        std::string(cTotalLength * (M + 1.) / nIter, '#') <<        // printing filled part
                        std::string(cTotalLength * (1. - (M + 1.) / nIter), '-') << // printing empty part
                        "] " << int((M + 1.) / nIter * 100.0) << "%\r";             // printing percentage

        // Updating regression coefficients and hyperparameters

        if (q > 0)
            PSBC::updateRP_clinical_cpp(p, q, x, ind_r, ind_d, ind_r_d, J, beta_prop_me, hyperPar_beta_prop_sd, xbeta, beta, h, sd_be, sampleRPc_accept);

        if (rw)
        {
            PSBC::updateRP_genomic_rw_cpp(p, x, ind_r, ind_d, ind_r_d, J, beta_prop_me, hyperPar_beta_prop_sd, xbeta, beta, h, sd_be, sampleRPg_accept);
        }
        else
        {
            PSBC::updateRP_genomic_cpp(p, x, ind_r, ind_d, ind_r_d, J, xbeta, beta, h, sd_be, sampleRPg_accept);
        }

        //        if( q > 0 )
        //        {
        for (unsigned int i = 0; i < K; ++i)
            be_normSq(i) = arma::accu(beta.elem(arma::find(groupInd == groupNo(i))) % beta.elem(arma::find(groupInd == groupNo(i))));
        //        } else {
        //            be_normSq = beta % beta;
        //        }

        h = PSBC::updateBH_cpp(ind_r_d, hPriorSh, d, hyperPar_c0, J, xbeta);

        tauSq = PSBC::updateTau_GL_cpp(hyperPar_lambdaSq, hyperPar_sigmaSq, be_normSq);

        hyperPar_sigmaSq = PSBC::updateSigma_GL_cpp(p, be_normSq, tauSq);
        hyperPar_lambdaSq = PSBC::updateLambda_GL_cpp(p, K, hyperPar_r, hyperPar_delta, tauSq);

        //        if( q > 0 )
        //        {
        for (unsigned int i = 0; i < K; ++i)
            tauSq_exp.elem(arma::find(groupInd == groupNo(i))).fill(tauSq(i));

        sd_bePart2.fill(hyperPar_beta_clin_sd);
        sd_be = arma::join_cols(sqrt(hyperPar_sigmaSq * tauSq_exp), sd_bePart2);
        //        } else {
        //            tauSq_exp = tauSq;
        //            sd_be = sqrt(sigmaSq * tauSq_exp);
        //        }

        // Save all results
        if (M % thin == 0)
        {
            ++j;
            beta_p.row(j) = beta.t();
            h_p.row(j) = h.t();
            tauSq_p.row(j) = tauSq.t();
            sigmaSq_p(j) = hyperPar_sigmaSq;
            lambdaSq_p(j) = hyperPar_lambdaSq;
        }

        // adaptive jumping rule
        if (j > 20)
            for (unsigned int jj = 0; jj < p + q; ++jj)
                if (beta(jj) == beta_p(j - 20, jj))
                    beta_prop_me(jj) = beta(jj);
    }

    arma::uvec accept_rate = join_cols(sampleRPg_accept, sampleRPc_accept); // / (double)(nIter);

    // Exit
    Rcpp::Rcout << "\nDONE, exiting! \n";

    return Rcpp::List::create(Rcpp::Named("beta.p") = beta_p,
                              Rcpp::Named("h.p") = h_p,
                              Rcpp::Named("tauSq.p") = tauSq_p,
                              Rcpp::Named("sigmaSq.p") = sigmaSq_p,
                              Rcpp::Named("lambdaSq.p") = lambdaSq_p,
                              Rcpp::Named("accept.rate") = accept_rate);
}
