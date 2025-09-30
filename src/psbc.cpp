#include "psbc.h"

PSBC::PSBC()
{
}

PSBC::~PSBC()
{
}

// TODO: maybe no need the following func
// multiply (element-wise) a matrix to a expanded vector
arma::mat PSBC::matProdVec(
    const arma::mat& x, 
    const arma::vec& y)
{
    // arma::mat mat_y = arma::zeros<arma::mat>(y.n_elem, x.n_cols);
    // mat_y.each_col() = y;
    // arma::mat spanMat = x % mat_y; // elementwise product
    // return spanMat;
    return x.each_col() % y;
}

// compute "arma::sum( matProdVec( ind_r_d_, exp_xbeta ).t(), 1 );"
arma::vec PSBC::sumMatProdVec(
    const arma::mat& x, 
    const arma::vec& y)
{
    // arma::vec spanVec = arma::zeros(x.n_cols);
    // for (unsigned int i = 0; i < x.n_cols; ++i)
    //     spanVec(i) = arma::dot(x.col(i), y);
    // return spanVec;
    return x.t() * y;
}

// TODO: test if passing (const) addresses are more efficient
// set a finite partition of the time axis to define the indicator matrices for risk sets and failure sets
// in order to calculate the increment in the cumulative baseline hazard in each interval
// also in order to construct the grouped data likelihood
void PSBC::settingInterval_cpp(
    const arma::vec& y, 
    const arma::uvec& delta_, 
    const arma::vec& s_, 
    const unsigned int J_,
    arma::mat &ind_d_,
    arma::mat &ind_r_,
    arma::mat &ind_r_d_,
    arma::vec &d_)
{
    // ind_d_ = ind_r_ = arma::zeros<arma::mat>(y.n_elem, J_);

    arma::uvec case0yleq;
    arma::uvec case0ygeq;
    arma::uvec case1yleq;
    arma::uvec case1ygeq;

    double smax = max(s_);
    case0yleq = arma::find(delta_ == 0 && y <= smax);
    case0ygeq = arma::find(delta_ == 0 && y > smax);
    case1yleq = arma::find(delta_ == 1 && y <= smax);
    case1ygeq = arma::find(delta_ == 1 && y > smax);

    int cen_j;
    for (unsigned int i = 0; i < case1yleq.n_elem; ++i)
    {
        cen_j = min(arma::find(s_ >= y(case1yleq(i))));
        ind_d_(case1yleq(i), cen_j) = 1.;
        ind_r_.submat(case1yleq(i), 0, case1yleq(i), cen_j).fill(1.);
        // Rcout << cen_j << ",";
    }

    for (unsigned int i = 0; i < case0yleq.n_elem; ++i)
    {
        cen_j = min(arma::find(s_ >= y(case0yleq(i))));
        ind_r_.submat(case0yleq(i), 0, case0yleq(i), cen_j).fill(1.);
    }

    if (case0ygeq.n_elem + case1ygeq.n_elem > 0)
    {
        arma::uvec union_case01ygeq;
        union_case01ygeq = arma::unique(arma::join_cols(case0ygeq, case1ygeq));
        ind_r_.rows(union_case01ygeq).fill(1.);
    }

    ind_r_d_ = ind_r_ - ind_d_;
    d_ = arma::sum(ind_d_.t(), 1);
}


// update cumulative baseline harzard
// update the increment h_j in the cumulative baseline hazard in each interval
arma::vec PSBC::updateBH_cpp(
    const arma::mat& ind_r_d_, 
    const arma::vec& hPriorSh_, 
    const arma::vec& d_, 
    double c0_, 
    const unsigned int J_, 
    const arma::vec& xbeta_)
{
    // arma::mat exp_xbeta_mat = matProdVec( ind_r_d_, arma::exp( xbeta_ ) );
    // arma::vec h_rate = c0_ + arma::sum( exp_xbeta_mat.t(), 1 );
    /*arma::vec h_rate( J_, arma::fill::value(c0_) ); // this is a bit faster than above
    xbeta_ = arma::exp( xbeta_ );
    for( unsigned int i=0; i<J_; ++i )
     //h_rate(i) += arma::accu( ind_r_d_.col(i) % arma::exp( xbeta_ ) );
        h_rate(i) += arma::dot( ind_r_d_.col(i), xbeta_ );*/
    arma::vec h_rate = c0_ + sumMatProdVec(ind_r_d_, arma::exp( xbeta_ ));

    // arma::vec shape = hPriorSh_ + d_;
    arma::vec h_ = arma::zeros<arma::vec>(J_);
    for (unsigned int j = 0; j < J_; ++j)
        // double h_rate = c0_ + arma::dot( ind_r_d_.col(j), arma::exp( xbeta_ ) );
        h_(j) = arma::randg(arma::distr_param(hPriorSh_(j) + d_(j), 1. / h_rate(j)));
    // h_(j) = R::rgamma( shape(j), 1. / h_rate(j) );

    return h_;
}

// sample values from an inverse-Gaussian distribution
arma::vec PSBC::rinvgauss(const arma::vec& a, double b)
{
    unsigned int n = a.n_elem;
    arma::vec pars(n);
    // Pre-generate all random numbers for better performance
    arma::vec z = arma::randn(n);  // z_i = R::rnorm(0., 1.);
    arma::vec u = arma::randu(n);  // u_i = R::runif(0., 1.);

    for (unsigned int i = 0; i < n; ++i)
    {
        double a_i = a(i); 
        double z_i = z(i);
        double u_i = u(i); 
        double y = z_i * z_i;
        double sqrt_term = std::sqrt(4.0 * a_i * b * y + a_i * a_i * y * y);
        double x = a_i + 0.5 * a_i * a_i * y / b - 0.5 * (a_i/ b) * sqrt_term;
        
        if (u_i <= a_i / (a_i + x))
        {
            pars(i) = x;
        }
        else
        {
            pars(i) = a_i * a_i / x;
        }
    }
    return pars;
}

// update hyperparameter tau (variance shrinkage of coefficients) sampled from the full conditional inverse-Gaussian distribution
arma::vec PSBC::updateTau_GL_cpp(
    double lambdaSq_, 
    double sigmaSq_, 
    const arma::vec& be_normSq_)
{
    arma::vec nu = arma::ones<arma::vec>(be_normSq_.n_elem);
    if (arma::any(be_normSq_ != 0))
    {
        nu = sqrt(lambdaSq_ * sigmaSq_ / be_normSq_);
        if (nu.has_inf())
            nu.elem(arma::find_nonfinite(nu)).fill(max(nu.elem(arma::find_finite(nu))) + 10.);
    }
    else
    {
        nu.fill(10.);
    }

    arma::vec tauSq = 1. / rinvgauss(nu, lambdaSq_);

    return tauSq;
}

// update variance parameter sigma_square sampled from the full conditional inverse-gamma distribution
double PSBC::updateSigma_GL_cpp(
    const unsigned int p, 
    const arma::vec& be_normSq_, 
    const arma::vec& tauSq_)
{
    double rate_sig = 0.5 * arma::accu(be_normSq_ / tauSq_);
    // if( rate_sig == 0. ) rate_sig = 0.0001;
    // double sigmaSq = 1. / R::rgamma( p / 2., 1. / rate_sig );
    double sigmaSq = 1. / arma::randg(arma::distr_param(p / 2., 1. / rate_sig));

    return sigmaSq;
}

// update hyperparameter lambda (variance shrinkage of tau) sampled from the full conditional gamma distribution
double PSBC::updateLambda_GL_cpp(
    const unsigned int p, 
    const unsigned int K, 
    double r, 
    double delta, 
    const arma::vec& tauSq_)
{
    double sumTauSq = arma::accu(tauSq_);
    double shape = (p + K) / 2. + r;
    double rate_lam = sumTauSq / 2. + delta;
    // double lambdaSq = R::rgamma( shape, 1. / rate );
    double lambdaSq = arma::randg(arma::distr_param(shape, 1. / rate_lam));

    return lambdaSq;
}

// update coefficients of clinical variables via a rw MH sampler
void PSBC::updateRP_clinical_cpp(
    const unsigned int p, 
    const unsigned int q, 
    const arma::mat& x_, 
    const arma::mat &ind_r_, 
    const arma::mat &ind_d_, 
    const arma::mat &ind_r_d_, 
    const unsigned int J_, 
    const arma::vec beta_prop_me_, 
    double beta_prop_sd, 
    arma::vec &xbeta_, 
    arma::vec &be_, 
    const arma::vec &h_, 
    const arma::vec &sd_be_, 
    arma::uvec &sampleRPc_accept_)
{
    // select parameters to be updated; use p+j for clinical
    arma::uvec updatej = arma::randperm(q);

    // declare relevant variables
    // double loglh_ini = 0.;
    // double loglh_prop = 0.;
    // double logprior_prop = 0.;
    // double logprior_ini = 0.;
    // double logprop_prop = 0.;
    // double logprop_ini = 0.;
    // double logR = 0.;

    // arma::vec be_prop;
    // arma::vec exp_xbeta;
    // arma::mat h_exp_xbeta_mat, h_exp_xbeta_prop_mat;
    // arma::vec first_sum, second_sum;
    // arma::vec first_sum_prop, second_sum_prop;
    // arma::vec xbeta_prop, exp_xbeta_prop;

    unsigned int j = 0;
    for (unsigned int j_id = 0; j_id < q; ++j_id)
    {
        j = updatej(j_id);
        // be_prop = be_;
        xbeta_.elem(arma::find(xbeta_ > 700)).fill(700.);
        arma::vec exp_xbeta = arma::exp(xbeta_);
        // first_sum = arma::sum( matProdVec(ind_r_d_, exp_xbeta).t(), 1 );
        arma::vec first_sum = sumMatProdVec(ind_r_d_, exp_xbeta);
        arma::mat h_exp_xbeta_mat = -arma::kron(exp_xbeta, h_.t());
        h_exp_xbeta_mat.elem(arma::find(h_exp_xbeta_mat > -1.0e-7)).fill(-1.0e-7);
        h_exp_xbeta_mat = arma::log(1.0 - arma::exp(h_exp_xbeta_mat));
        arma::vec second_sum = arma::sum((h_exp_xbeta_mat % ind_d_).t(), 1);
        double loglh_ini = arma::accu(-h_ % first_sum + second_sum);

        // clinical version:
        // be_prop( p + j ) = R::rnorm( beta_prop_me_(p+j), beta_prop_sd );
        double be_prop = arma::randn(arma::distr_param(beta_prop_me_(p + j), beta_prop_sd));
        arma::vec xbeta_prop = xbeta_ - x_.col(p + j) * be_(p + j) + x_.col(p + j) * be_prop;
        xbeta_prop.elem(arma::find(xbeta_prop > 700.)).fill(700.);
        arma::vec exp_xbeta_prop = arma::exp(xbeta_prop);
        // first_sum_prop = arma::sum( matProdVec( ind_r_d_, exp_xbeta_prop ).t(), 1);
        arma::vec first_sum_prop = sumMatProdVec(ind_r_d_, exp_xbeta_prop);

        arma::mat h_exp_xbeta_prop_mat = -arma::kron(exp_xbeta_prop, h_.t());
        h_exp_xbeta_prop_mat.elem(arma::find(h_exp_xbeta_prop_mat > -1.0e-7)).fill(-1.0e-7);
        h_exp_xbeta_prop_mat = arma::log(1.0 - arma::exp(h_exp_xbeta_prop_mat));
        arma::vec second_sum_prop = arma::sum((h_exp_xbeta_prop_mat % ind_d_).t(), 1);

        double loglh_prop = arma::accu(-h_ % first_sum_prop + second_sum_prop);
        /*logprior_prop = R::dnorm( be_prop(p+j), 0.0, sd_be_(p+j), true);
        logprior_ini = R::dnorm( be_(p+j), 0.0, sd_be_(p+j), true);
        logprop_prop = R::dnorm( be_prop(p+j), beta_prop_me_(p+j), beta_prop_sd, true);
        logprop_ini = R::dnorm( be_(p+j), beta_prop_me_(p+j), beta_prop_sd, true);*/
        double logprior_prop = arma::log_normpdf(be_prop, 0.0, sd_be_(p + j));
        double logprior_ini = arma::log_normpdf(be_(p + j), 0.0, sd_be_(p + j));
        double logprop_prop = arma::log_normpdf(be_prop, beta_prop_me_(p + j), beta_prop_sd);
        double logprop_ini = arma::log_normpdf(be_(p + j), beta_prop_me_(p + j), beta_prop_sd);
        double logR = loglh_prop - loglh_ini + logprior_prop - logprior_ini + logprop_ini - logprop_prop;

        // if( log( R::runif(0., 1.) ) < logR )
        if (log(arma::randu()) < logR)
        {
            be_(p + j) = be_prop;
            xbeta_ = xbeta_prop;
            sampleRPc_accept_(j)++;
        }
    }
}

// update coefficients of genomic variables via a MH sampler
void PSBC::updateRP_genomic_cpp(
    const unsigned int p, 
    const arma::mat& x_, 
    const arma::mat &ind_r_, 
    const arma::mat &ind_d_, 
    const arma::mat &ind_r_d_, 
    const unsigned int J_, 
    arma::vec &xbeta_, 
    arma::vec &be_, 
    const arma::vec &h_, 
    const arma::vec &sd_be_, 
    arma::uvec &sampleRPg_accept_)
{
    arma::uvec updatej = arma::randperm(p);

    unsigned int j = 0;
    for (unsigned int j_id = 0; j_id < p; ++j_id)
    {
        j = updatej(j_id);
        xbeta_.elem(arma::find(xbeta_ > 700)).fill(700.);
        arma::vec exp_xbeta = arma::exp(xbeta_);
        arma::vec x_exp_xbeta = x_.col(j) % exp_xbeta;
        // D1_1st = - h_ % arma::sum( matProdVec( ind_r_d_, x_exp_xbeta ).t(), 1 );
        arma::vec D1_1st = -h_ % sumMatProdVec(ind_r_d_, x_exp_xbeta);

        arma::mat h_exp_xbeta_mat = -arma::kron(exp_xbeta, h_.t());
        h_exp_xbeta_mat.elem(arma::find(h_exp_xbeta_mat > -1.0e-7)).fill(-1.0e-7);
        arma::mat exp_h_exp_xbeta_mat = arma::exp(h_exp_xbeta_mat);
        arma::mat D1_2nd_den = 1. - exp_h_exp_xbeta_mat;
        arma::mat D1_2nd_num = matProdVec(exp_h_exp_xbeta_mat, x_exp_xbeta);
        arma::vec D1_2nd = h_ % arma::sum((D1_2nd_num / D1_2nd_den % ind_d_).t(), 1);
        double D1 = arma::sum(D1_1st + D1_2nd) - 1. / sd_be_(j) / sd_be_(j) * be_(j);

        arma::vec x_sq_exp_xbeta = x_.col(j) % x_.col(j) % exp_xbeta;
        // D2_1st = - h_ % arma::sum( matProdVec( ind_r_d_, x_sq_exp_xbeta ).t(), 1 );
        arma::vec D2_1st = -h_ % sumMatProdVec(ind_r_d_, x_sq_exp_xbeta);
        arma::mat D2_2nd_den = D1_2nd_den % D1_2nd_den;
        arma::mat D2_2nd_num = matProdVec(exp_h_exp_xbeta_mat, x_sq_exp_xbeta) % (1. - exp_h_exp_xbeta_mat + h_exp_xbeta_mat);
        arma::vec D2_2nd = h_ % arma::sum((D2_2nd_num / D2_2nd_den % ind_d_).t(), 1);
        double D2 = arma::accu(D2_1st + D2_2nd) - 1. / sd_be_(j) / sd_be_(j);

        double be_prop_me = be_(j) - D1 / D2;
        double be_prop_sd = 2.4 / sqrt(-D2);
        // be_prop = be_;

        // genomic version:
        // be_prop = R::rnorm( be_prop_me, be_prop_sd );
        double be_prop = arma::randn(arma::distr_param(be_prop_me, be_prop_sd));
        arma::vec xbeta_prop = xbeta_ - x_.col(j) * be_(j) + x_.col(j) * be_prop;
        xbeta_prop.elem(arma::find(xbeta_prop > 700)).fill(700.);
        arma::vec exp_xbeta_prop = arma::exp(xbeta_prop);
        arma::vec x_exp_xbeta_prop = x_.col(j) % exp_xbeta_prop;
        // D1_1st_prop = - h_ % arma::sum( matProdVec( ind_r_d_, x_exp_xbeta_prop ).t(), 1 );
        arma::vec D1_1st_prop = -h_ % sumMatProdVec(ind_r_d_, x_exp_xbeta_prop);

        arma::mat h_exp_xbeta_prop_mat = -arma::kron(exp_xbeta_prop, h_.t());
        h_exp_xbeta_prop_mat.elem(arma::find(h_exp_xbeta_prop_mat > -1.0e-7)).fill(-1.0e-7);
        arma::mat exp_h_exp_xbeta_prop_mat = arma::exp(h_exp_xbeta_prop_mat);
        arma::mat D1_2nd_den_prop = 1. - exp_h_exp_xbeta_prop_mat;
        arma::mat D1_2nd_num_prop = matProdVec(exp_h_exp_xbeta_prop_mat, x_exp_xbeta_prop);
        arma::vec D1_2nd_prop = h_ % arma::sum((D1_2nd_num_prop / D1_2nd_den_prop % ind_d_).t(), 1);
        double D1_prop = arma::accu(D1_1st_prop + D1_2nd_prop) - 1. / sd_be_(j) / sd_be_(j) * be_prop;

        arma::vec x_sq_exp_xbeta_prop = x_.col(j) % x_.col(j) % exp_xbeta_prop;
        // D2_1st_prop = -h_ % arma::sum( matProdVec( ind_r_d_, x_sq_exp_xbeta_prop ).t(), 1);
        arma::vec D2_1st_prop = -h_ % sumMatProdVec(ind_r_d_, x_sq_exp_xbeta_prop);
        arma::mat D2_2nd_den_prop = D1_2nd_den_prop % D1_2nd_den_prop;
        arma::mat D2_2nd_num_prop = matProdVec(exp_h_exp_xbeta_prop_mat, x_sq_exp_xbeta_prop) % (1. - exp_h_exp_xbeta_prop_mat + h_exp_xbeta_prop_mat);
        arma::vec D2_2nd_prop = h_ % arma::sum((D2_2nd_num_prop / D2_2nd_den_prop, ind_d_).t(), 1);
        double D2_prop = arma::accu(D2_1st_prop + D2_2nd_prop) - 1. / sd_be_(j) / sd_be_(j);
        double be_prop_me_ini = be_prop - D1_prop / D2_prop;
        double be_prop_sd_ini = 2.4 / sqrt(-D2_prop);

        // first_sum = arma::sum( matProdVec( ind_r_d_, exp_xbeta ).t(), 1 );
        arma::vec first_sum = sumMatProdVec(ind_r_d_, exp_xbeta);
        arma::vec second_sum = arma::sum((arma::log(D1_2nd_den) % ind_d_).t(), 1);

        double loglh_ini = arma::accu(-h_ % first_sum + second_sum);
        // first_sum_prop = arma::sum( matProdVec( ind_r_d_, exp_xbeta_prop ).t(), 1) ;
        arma::vec first_sum_prop = sumMatProdVec(ind_r_d_, exp_xbeta_prop);
        arma::vec second_sum_prop = arma::sum((arma::log(D1_2nd_den_prop) % ind_d_).t(), 1);
        double loglh_prop = arma::accu(-h_ % first_sum_prop + second_sum_prop);

        /*logprior_prop = R::dnorm( be_prop, 0.0, sd_be_(j), true);
        logprior_ini = R::dnorm( be_(j), 0.0, sd_be_(j), true);
        logprop_prop = R::dnorm( be_prop, be_prop_me_ini, be_prop_sd_ini, true);
        logprop_ini = R::dnorm( be_(j), be_prop_me, be_prop_sd, true);*/
        double logprior_prop = arma::log_normpdf(be_prop, 0.0, sd_be_(j));
        double logprior_ini = arma::log_normpdf(be_(j), 0.0, sd_be_(j));
        double logprop_prop = arma::log_normpdf(be_prop, be_prop_me_ini, be_prop_sd_ini);
        double logprop_ini = arma::log_normpdf(be_(j), be_prop_me, be_prop_sd);
        double logR = loglh_prop - loglh_ini + logprior_prop - logprior_ini + logprop_ini - logprop_prop;

        // if( log( R::runif(0., 1.) ) < logR )
        if (log(arma::randu()) < logR)
        {
            be_(j) = be_prop;
            xbeta_ = xbeta_prop;
            sampleRPg_accept_(j)++;
        }
    }
}

// update coefficients of genomic variables via a rw MH sampler, almost the same as updateRP_clinical_cpp()
void PSBC::updateRP_genomic_rw_cpp(
    const unsigned int p, 
    const arma::mat& x_, 
    const arma::mat &ind_r_, 
    const arma::mat &ind_d_, 
    const arma::mat &ind_r_d_, 
    const unsigned int J_, 
    const arma::vec beta_prop_me_, 
    double beta_prop_sd, 
    arma::vec &xbeta_, 
    arma::vec &be_, 
    const arma::vec &h_, 
    const arma::vec &sd_be_, 
    arma::uvec &sampleRPg_accept_)
{
    // select parameters to be updated; use p+j for clinical
    arma::uvec updatej = arma::randperm(p);

    unsigned int j = 0;
    for (unsigned int j_id = 0; j_id < p; ++j_id)
    {
        j = updatej(j_id);
        // be_prop = be_;
        xbeta_.elem(arma::find(xbeta_ > 700)).fill(700.);
        arma::vec exp_xbeta = arma::exp(xbeta_);
        // first_sum = arma::sum( matProdVec(ind_r_d_, exp_xbeta).t(), 1 );
        arma::vec first_sum = sumMatProdVec(ind_r_d_, exp_xbeta);
        arma::vec h_exp_xbeta_mat = -arma::kron(exp_xbeta, h_.t());
        h_exp_xbeta_mat.elem(arma::find(h_exp_xbeta_mat > -1.0e-7)).fill(-1.0e-7);
        h_exp_xbeta_mat = arma::log(1.0 - arma::exp(h_exp_xbeta_mat));
        arma::vec second_sum = arma::sum((h_exp_xbeta_mat % ind_d_).t(), 1);
        double loglh_ini = arma::accu(-h_ % first_sum + second_sum);

        // clinical version:
        // be_prop( j ) = R::rnorm( beta_prop_me_(j), beta_prop_sd );
        double be_prop = arma::randn(arma::distr_param(beta_prop_me_(j), beta_prop_sd));
        arma::vec xbeta_prop = xbeta_ - x_.col(j) * be_(j) + x_.col(j) * be_prop;
        xbeta_prop.elem(arma::find(xbeta_prop > 700)).fill(700.);
        arma::vec exp_xbeta_prop = arma::exp(xbeta_prop);
        // first_sum_prop = arma::sum( matProdVec( ind_r_d_, exp_xbeta_prop ).t(), 1);
        arma::vec first_sum_prop = sumMatProdVec(ind_r_d_, exp_xbeta_prop);

        arma::vec h_exp_xbeta_prop_mat = -arma::kron(exp_xbeta_prop, h_.t());
        h_exp_xbeta_prop_mat.elem(arma::find(h_exp_xbeta_prop_mat > -1.0e-7)).fill(-1.0e-7);
        h_exp_xbeta_prop_mat = 1.0 - arma::exp(h_exp_xbeta_prop_mat);
        arma::vec second_sum_prop = arma::sum((arma::log(h_exp_xbeta_prop_mat) % ind_d_).t(), 1);

        double loglh_prop = arma::accu(-h_ % first_sum_prop + second_sum_prop);
        /*logprior_prop = R::dnorm( be_prop, 0.0, sd_be_(j), true);
        logprior_ini = R::dnorm( be_(j), 0.0, sd_be_(j), true);
        logprop_prop = R::dnorm( be_prop, beta_prop_me_(j), beta_prop_sd, true);
        logprop_ini = R::dnorm( be_(j), beta_prop_me_(j), beta_prop_sd, true);*/
        double logprior_prop = arma::log_normpdf(be_prop, 0.0, sd_be_(j));
        double logprior_ini = arma::log_normpdf(be_(j), 0.0, sd_be_(j));
        double logprop_prop = arma::log_normpdf(be_prop, beta_prop_me_(j), beta_prop_sd);
        double logprop_ini = arma::log_normpdf(be_(j), beta_prop_me_(j), beta_prop_sd);
        double logR = loglh_prop - loglh_ini + logprior_prop - logprior_ini + logprop_ini - logprop_prop;

        // if( log( R::runif(0., 1.) ) < logR )
        if (log(arma::randu()) < logR)
        {
            be_(j) = be_prop;
            xbeta_ = xbeta_prop;
            sampleRPg_accept_(j)++;
        }
    }
}
