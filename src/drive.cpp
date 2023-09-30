#include "drive.h"

#ifndef CCODE
    using Rcpp::Rcout;
    using Rcpp::Rcerr;
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

// convert a vector into a matrix by copying columns
arma::mat vecToMat( const arma::vec x, const int ncol )
{
    arma::mat spanMat = arma::zeros<arma::mat>( x.n_elem, ncol );
    spanMat.each_col() = x;
    return spanMat;
}

// multiply (element-wise) a matrix to a reshaped vector (via function vecToMat())
arma::mat matProdVec( const arma::mat x, const arma::vec& y )
{
    arma::mat mat_y = arma::zeros<arma::mat>( y.n_elem, x.n_cols );
    mat_y.each_col() = y;
    arma::mat spanMat =  x % mat_y; // elementwise product
    return spanMat;
}

// set a finite partition of the time axis to define the indicator matrices for risk sets and failure sets
// in order to calculate the increment in the cumulative baseline hazard in each interval
// also in order to construct the grouped data likelihood
void settingInterval_cpp( const arma::vec y, const arma::vec delta_, const arma::vec s_, unsigned int J_, arma::mat& ind_d_, arma::mat& ind_r_, arma::mat& ind_r_d_, arma::vec& d_ )
{
    ind_d_ = ind_r_ = arma::zeros<arma::mat>(y.n_elem, J_);
    
    arma::uvec case0yleq;
    arma::uvec case0ygeq;
    arma::uvec case1yleq;
    arma::uvec case1ygeq;
    
    double smax = max(s_);
    case0yleq = arma::find( delta_ == 0. && y <= smax );
    case0ygeq = arma::find( delta_ == 0. && y > smax );
    case1yleq = arma::find( delta_ == 1. && y <= smax );
    case1ygeq = arma::find( delta_ == 1. && y > smax );
  
    int cen_j;
    for( unsigned int i = 0; i < case1yleq.n_elem; ++i )
    {
        cen_j = min( arma::find( s_ >= y(case1yleq(i)) ) );
        ind_d_( case1yleq(i), cen_j ) = 1.;
        ind_r_.submat( case1yleq(i), 0, case1yleq(i), cen_j ).fill( 1. );
        //Rcout << cen_j << ",";
    }
    
    for( unsigned int i = 0; i < case0yleq.n_elem; ++i )
    {
        cen_j = min( arma::find( s_ >= y(case0yleq(i)) ) );
        ind_r_.submat( case0yleq(i), 0, case0yleq(i), cen_j ).fill( 1. );
    }
    
    if( case0ygeq.n_elem + case1ygeq.n_elem > 0 )
    {
        arma::uvec union_case01ygeq;
        union_case01ygeq = arma::unique( arma::join_cols(case0ygeq, case1ygeq) );
        ind_r_.rows( union_case01ygeq ).fill( 1. );
    }
    
    ind_r_d_ = ind_r_ - ind_d_;
    d_ = arma::sum( ind_d_.t(), 1 );
}

// update cumulative baseline harzard
// update the increment h_j in the cumulative baseline hazard in each interval
arma::vec updateBH_cpp( arma::mat& ind_r_d_, arma::vec hPriorSh_, arma::vec& d_, double c0_, unsigned int J_, arma::vec xbeta_ )
{
    arma::mat exp_xbeta_mat = matProdVec( ind_r_d_, arma::exp( xbeta_ ) );
    arma::vec h_rate = c0_ + arma::sum( exp_xbeta_mat.t(), 1 );
    arma::vec shape = hPriorSh_ + d_;
 
    arma::vec h_ = arma::zeros<arma::vec>( J_ );
    for( unsigned int j=0; j<J_; ++j )
        h_(j) = arma::randg( arma::distr_param( shape(j), 1. / h_rate(j) ) );
//        if( shape(j) != 0. )
//        {
//            //h_(j) = R::rgamma( shape(j), 1. / h_rate(j) );
//            h_(j) = arma::randg( arma::distr_param( shape(j), 1. / h_rate(j) ) );
//        } else {
//            h_(j) = 0; // since R::rgamma( 0, 1./h_rate(j) ) = 0
//        }
    
    return h_;
}

// sample values from an inverse-Gaussian distribution
arma::vec rinvgauss( arma::vec a, double b )
{
    int n = a.n_elem;
    arma::vec pars = arma::zeros<arma::vec>(n);
    double z, y, x, u;

    for( unsigned int i=0; i<n; ++i )
    {
        //z = R::rnorm(0., 1.);
        z = arma::randn();
        y = z*z;
        x = a(i) + 0.5*a(i)*a(i)*y/b - 0.5 * (a(i)/b) * sqrt( 4.*a(i)*b*y+a(i)*a(i)*y*y );
        //u = R::runif(0., 1.);
        u = arma::randu();
        if( u <= a(i) / ( a(i)+x ) )
        {
            pars(i) = x;
        } else {
            pars(i) = a(i) * a(i) / x;
        }
    }
    return pars;
}

// update hyperparameter tau (variance shrinkage of coefficients) sampled from the full conditional inverse-Gaussian distribution
arma::vec updateTau_GL_cpp( double lambdaSq_, double sigmaSq_, arma::vec be_normSq_ )
{
    arma::vec nu = arma::ones<arma::vec>( be_normSq_.n_elem );
    if( arma::any( be_normSq_ != 0 ) )
    {
        nu = sqrt( lambdaSq_ * sigmaSq_ / be_normSq_ );
        if( nu.has_inf() )
            nu.elem( arma::find_nonfinite(nu)).fill( max( nu.elem(arma::find_finite(nu)) ) + 10. );
        
    } else {
        nu.fill( 10. );
    }
    
    arma::vec tauSq = rinvgauss( nu, lambdaSq_ );
    
    return tauSq;
}


// update variance parameter sigma_square sampled from the full conditional inverse-gamma distribution
double updateSigma_GL_cpp( int p, arma::vec be_normSq_, arma::vec tauSq_ )
{
    double rate_sig = 0.5 * arma::accu( be_normSq_ / tauSq_ );
    //if( rate_sig == 0. ) rate_sig = 0.0001;
    //double sigmaSq = 1. / R::rgamma( p / 2., 1. / rate_sig );
    double sigmaSq = 1. / arma::randg( arma::distr_param( p / 2., 1. / rate_sig ) );
    
    return sigmaSq;
}

// update hyperparameter lambda (variance shrinkage of tau) sampled from the full conditional gamma distribution
double updateLambda_GL_cpp( int p, unsigned int K, double r, double delta, arma::vec tauSq_ )
{
    double sumTauSq = arma::accu( tauSq_ );
    double shape = (p + K) / 2. + r;
    double rate_lam = sumTauSq / 2. + delta;
    //double lambdaSq = R::rgamma( shape, 1. / rate );
    double lambdaSq = arma::randg( arma::distr_param( shape, 1. / rate_lam ) );
    
    return lambdaSq;
}

// update coefficients of clinical variables via a rw MH sampler
void updateRP_clinical_cpp( int p, int q, const arma::mat x_, arma::mat& ind_r_, arma::mat& ind_d_, arma::mat& ind_r_d_, unsigned int J_, arma::vec beta_prop_me_, double beta_prop_sd, arma::vec& xbeta_, arma::vec& be_, arma::vec& h_, arma::vec sd_be_, arma::uvec& sampleRPc_accept_ )
{
    // select parameters to be updated; use p+j for clinical
    arma::uvec updatej = arma::randperm( q );
    
    arma::vec be_prop = be_;
    arma::vec exp_xbeta = arma::exp( xbeta_ );
    arma::vec xbeta_prop;
    arma::vec exp_xbeta_prop;
    
    arma::mat h_exp_xbeta_mat, h_exp_xbeta_prop_mat;
    
    arma::vec first_sum, second_sum;
    arma::vec first_sum_prop, second_sum_prop;
    double loglh_ini, loglh_prop, logprior_prop, logprior_ini, logprop_prop, logR, logprop_ini;
    
    int j = 0;
    for( int j_id=0; j_id<q; ++j_id )
    {
        j = updatej(j_id);
        be_prop = be_;
        xbeta_.elem( arma::find(xbeta_ > 700) ).fill( 700. );
        exp_xbeta = arma::exp( xbeta_ );
        first_sum = arma::sum( matProdVec(ind_r_d_, exp_xbeta).t(), 1 );
        h_exp_xbeta_mat = - arma::kron( exp_xbeta, h_.t() );
        h_exp_xbeta_mat.elem( arma::find(h_exp_xbeta_mat > -1.0e-7) ).fill( -1.0e-7 );
        h_exp_xbeta_mat = arma::log( 1.0 - arma::exp(h_exp_xbeta_mat) );
        second_sum = arma::sum( ( h_exp_xbeta_mat % ind_d_ ).t(), 1 );
        loglh_ini = arma::accu( - h_ % first_sum + second_sum );
        
        // clinical version:
        //be_prop( p + j ) = R::rnorm( beta_prop_me_(p+j), beta_prop_sd );
        be_prop( p + j ) = arma::randn( arma::distr_param( beta_prop_me_(p+j), beta_prop_sd ) );
        xbeta_prop = xbeta_ - x_.col(p+j) * be_(p+j) + x_.col(p+j) * be_prop(p+j);
        xbeta_prop.elem( arma::find( xbeta_prop > 700.) ).fill( 700. );
        exp_xbeta_prop = arma::exp( xbeta_prop );
        first_sum_prop = arma::sum( matProdVec( ind_r_d_, exp_xbeta_prop ).t(), 1);
        
        h_exp_xbeta_prop_mat = - arma::kron( exp_xbeta_prop, h_.t() );
        h_exp_xbeta_prop_mat.elem( arma::find(h_exp_xbeta_prop_mat > -1.0e-7) ).fill( -1.0e-7 );
        h_exp_xbeta_prop_mat = arma::log( 1.0 - arma::exp(h_exp_xbeta_prop_mat) );
        second_sum_prop = arma::sum( ( h_exp_xbeta_prop_mat % ind_d_ ).t(), 1 );
        
        loglh_prop = arma::accu( - h_ % first_sum_prop + second_sum_prop );
        /*logprior_prop = R::dnorm( be_prop(p+j), 0.0, sd_be_(p+j), true);
        logprior_ini = R::dnorm( be_(p+j), 0.0, sd_be_(p+j), true);
        logprop_prop = R::dnorm( be_prop(p+j), beta_prop_me_(p+j), beta_prop_sd, true);
        logprop_ini = R::dnorm( be_(p+j), beta_prop_me_(p+j), beta_prop_sd, true);*/
        logprior_prop = arma::log_normpdf( be_prop(p+j), 0.0, sd_be_(p+j) );
        logprior_ini = arma::log_normpdf( be_(p+j), 0.0, sd_be_(p+j) );
        logprop_prop = arma::log_normpdf( be_prop(p+j), beta_prop_me_(p+j), beta_prop_sd );
        logprop_ini = arma::log_normpdf( be_(p+j), beta_prop_me_(p+j), beta_prop_sd );
        logR = loglh_prop - loglh_ini + logprior_prop - logprior_ini + logprop_ini - logprop_prop;
        
        //if( log( R::runif(0., 1.) ) < logR )
        if( log( arma::randu() ) < logR )
        {
            be_(p+j) = be_prop(p+j);
            xbeta_ = xbeta_prop;
            sampleRPc_accept_(j)++;
        }
    }
}

// update coefficients of genomic variables via a MH sampler
void updateRP_genomic_cpp( int p, const arma::mat x_, arma::mat& ind_r_, arma::mat& ind_d_, arma::mat& ind_r_d_, unsigned int J_, arma::vec& xbeta_, arma::vec& be_, arma::vec& h_, arma::vec sd_be_, arma::uvec& sampleRPg_accept_ )
{
    arma::uvec updatej = arma::randperm(p);
    
    arma::vec be_prop = be_;
    arma::vec exp_xbeta = arma::exp( xbeta_ );
    arma::mat h_exp_xbeta_mat, h_exp_xbeta_prop_mat;
    arma::mat exp_h_exp_xbeta_mat, exp_h_exp_xbeta_prop_mat;
    arma::vec x_sq_exp_xbeta, first_sum, first_sum_prop, second_sum, second_sum_prop;
    arma::vec x_exp_xbeta, xbeta_prop, exp_xbeta_prop, x_exp_xbeta_prop, x_sq_exp_xbeta_prop;
    arma::vec D1_1st, D1_2nd, D1_1st_prop, D1_2nd_prop;
    arma::vec D2_1st, D2_2nd, D2_1st_prop, D2_2nd_prop;
    arma::mat D1_2nd_den, D1_2nd_num, D1_2nd_den_prop, D1_2nd_num_prop;
    arma::mat D2_2nd_num, D2_2nd_den, D2_2nd_den_prop, D2_2nd_num_prop;
    
    //Rcout << "updateRP_genomic_cpp...test1\n";
    double be_prop_me = 1.;
    double be_prop_sd = 1.;
    double be_prop_me_ini, be_prop_sd_ini, D1, D2, D1_prop, D2_prop, loglh_ini, loglh_prop, logprior_prop, logprior_ini, logprop_prop, logprop_ini, logR;
    int j = 0;
    for( int j_id=0; j_id<p; ++j_id )
    {
        j = updatej(j_id);
        xbeta_.elem( arma::find(xbeta_ > 700) ).fill(700.);
        exp_xbeta = arma::exp(xbeta_);
        x_exp_xbeta = x_.col(j) % exp_xbeta;
        D1_1st = - h_ % arma::sum( matProdVec( ind_r_d_, x_exp_xbeta ).t(), 1 );
        
        h_exp_xbeta_mat = - arma::kron( exp_xbeta, h_.t() );
        h_exp_xbeta_mat.elem( arma::find(h_exp_xbeta_mat > -1.0e-7) ).fill(-1.0e-7);
        exp_h_exp_xbeta_mat = arma::exp( h_exp_xbeta_mat );
        D1_2nd_den = 1. - exp_h_exp_xbeta_mat;
        D1_2nd_num = matProdVec( exp_h_exp_xbeta_mat, x_exp_xbeta );
        D1_2nd = h_ % arma::sum( ( D1_2nd_num / D1_2nd_den % ind_d_ ).t(), 1 );
        D1 = arma::sum( D1_1st + D1_2nd ) - 1. /  sd_be_(j) / sd_be_(j) * be_(j);
        
        x_sq_exp_xbeta = x_.col(j) % x_.col(j) % exp_xbeta;
        D2_1st = - h_ % arma::sum( matProdVec( ind_r_d_, x_sq_exp_xbeta ).t(), 1 );
        D2_2nd_den = D1_2nd_den % D1_2nd_den;
        D2_2nd_num = matProdVec( exp_h_exp_xbeta_mat, x_sq_exp_xbeta ) % ( 1. - exp_h_exp_xbeta_mat + h_exp_xbeta_mat );
        D2_2nd = h_ % arma::sum( ( D2_2nd_num / D2_2nd_den % ind_d_ ).t(), 1 );
        D2 = arma::accu(D2_1st + D2_2nd) - 1. / sd_be_(j) / sd_be_(j);
        
        be_prop_me = be_(j) - D1 / D2;
        be_prop_sd = 2.4 / sqrt(-D2);
        be_prop = be_;
        
        // genomic version:
        //be_prop(j) = R::rnorm( be_prop_me, be_prop_sd );
        be_prop( j ) = arma::randn( arma::distr_param( be_prop_me, be_prop_sd ) );
        xbeta_prop = xbeta_ - x_.col(j) * be_(j) + x_.col(j) * be_prop(j);
        xbeta_prop.elem( arma::find(xbeta_prop > 700) ).fill( 700. );
        exp_xbeta_prop = arma::exp( xbeta_prop );
        x_exp_xbeta_prop = x_.col(j) % exp_xbeta_prop;
        D1_1st_prop = - h_ % arma::sum( matProdVec( ind_r_d_, x_exp_xbeta_prop ).t(), 1 );
        
        h_exp_xbeta_prop_mat = - arma::kron( exp_xbeta_prop, h_.t() );
        h_exp_xbeta_prop_mat.elem( arma::find(h_exp_xbeta_prop_mat > -1.0e-7) ).fill( -1.0e-7 );
        exp_h_exp_xbeta_prop_mat = arma::exp(h_exp_xbeta_prop_mat);
        D1_2nd_den_prop = 1. - exp_h_exp_xbeta_prop_mat;
        D1_2nd_num_prop = matProdVec( exp_h_exp_xbeta_prop_mat, x_exp_xbeta_prop );
        D1_2nd_prop = h_ % arma::sum( ( D1_2nd_num_prop / D1_2nd_den_prop % ind_d_ ).t(), 1 );
        D1_prop = arma::accu( D1_1st_prop + D1_2nd_prop ) - 1. / sd_be_(j) / sd_be_(j) * be_prop(j);
        
        x_sq_exp_xbeta_prop = x_.col(j) % x_.col(j) % exp_xbeta_prop;
        D2_1st_prop = -h_ % arma::sum( matProdVec( ind_r_d_, x_sq_exp_xbeta_prop ).t(), 1);
        D2_2nd_den_prop = D1_2nd_den_prop % D1_2nd_den_prop;
        D2_2nd_num_prop = matProdVec( exp_h_exp_xbeta_prop_mat, x_sq_exp_xbeta_prop ) % (1. - exp_h_exp_xbeta_prop_mat + h_exp_xbeta_prop_mat );
        D2_2nd_prop = h_ % arma::sum( ( D2_2nd_num_prop / D2_2nd_den_prop, ind_d_ ).t(), 1 );
        D2_prop = arma::accu( D2_1st_prop + D2_2nd_prop) - 1. / sd_be_(j) / sd_be_(j);
        be_prop_me_ini = be_prop(j) - D1_prop / D2_prop;
        be_prop_sd_ini = 2.4 / sqrt(-D2_prop);
        
        first_sum = arma::sum( matProdVec( ind_r_d_, exp_xbeta ).t(), 1 );
        second_sum = arma::sum( ( arma::log(D1_2nd_den) % ind_d_ ).t(), 1 );
        
        loglh_ini = arma::accu( - h_ % first_sum + second_sum );
        first_sum_prop = arma::sum( matProdVec( ind_r_d_, exp_xbeta_prop ).t(), 1) ;
        second_sum_prop = arma::sum( ( arma::log(D1_2nd_den_prop) % ind_d_ ).t(), 1 );
        loglh_prop = arma::accu( - h_ % first_sum_prop + second_sum_prop );
        
        /*logprior_prop = R::dnorm( be_prop(j), 0.0, sd_be_(j), true);
        logprior_ini = R::dnorm( be_(j), 0.0, sd_be_(j), true);
        logprop_prop = R::dnorm( be_prop(j), be_prop_me_ini, be_prop_sd_ini, true);
        logprop_ini = R::dnorm( be_(j), be_prop_me, be_prop_sd, true);*/
        logprior_prop = arma::log_normpdf( be_prop(j), 0.0, sd_be_(j) );
        logprior_ini = arma::log_normpdf( be_(j), 0.0, sd_be_(j) );
        logprop_prop = arma::log_normpdf( be_prop(j), be_prop_me_ini, be_prop_sd_ini );
        logprop_ini = arma::log_normpdf( be_(j), be_prop_me, be_prop_sd );
        logR = loglh_prop - loglh_ini + logprior_prop - logprior_ini + logprop_ini - logprop_prop;
        
        //if( log( R::runif(0., 1.) ) < logR )
        if( log( arma::randu() ) < logR )
        {
            be_(j) = be_prop( j );
            xbeta_ = xbeta_prop;
            sampleRPg_accept_( j ) = sampleRPg_accept_( j ) + 1;
        }
    }
}

// update coefficients of genomic variables via a rw MH sampler, almost the same as updateRP_clinical_cpp()
void updateRP_genomic_rw_cpp( int p, const arma::mat x_, arma::mat& ind_r_, arma::mat& ind_d_, arma::mat& ind_r_d_, unsigned int J_, arma::vec beta_prop_me_, double beta_prop_sd, arma::vec& xbeta_, arma::vec& be_, arma::vec& h_, arma::vec sd_be_, arma::uvec& sampleRPg_accept_ )
{
    // select parameters to be updated; use p+j for clinical
    arma::uvec updatej = arma::randperm( p );
    
    arma::vec be_prop = be_;
    arma::vec exp_xbeta = arma::exp( xbeta_ );
    arma::vec xbeta_prop;
    arma::vec exp_xbeta_prop;
    arma::mat h_exp_xbeta_mat, h_exp_xbeta_prop_mat;
    
    //Rcout << "updateRP_clinical_cpp...test1\n";
    arma::vec first_sum, second_sum;
    arma::vec first_sum_prop, second_sum_prop;
    double loglh_ini, loglh_prop, logprior_prop, logprior_ini, logprop_prop, logR, logprop_ini;
    
    int j = 0;
    for( int j_id=0; j_id<p; ++j_id )
    {
        j = updatej(j_id);
        be_prop = be_;
        xbeta_.elem( arma::find(xbeta_ > 700) ).fill( 700. );
        exp_xbeta = arma::exp( xbeta_ );
        first_sum = arma::sum( matProdVec(ind_r_d_, exp_xbeta).t(), 1 );
        h_exp_xbeta_mat = - arma::kron( exp_xbeta, h_.t() );
        h_exp_xbeta_mat.elem( arma::find(h_exp_xbeta_mat > -1.0e-7) ).fill( -1.0e-7 );
        h_exp_xbeta_mat = arma::log( 1.0 - arma::exp(h_exp_xbeta_mat) );
        second_sum = arma::sum( ( h_exp_xbeta_mat % ind_d_ ).t(), 1 );
        loglh_ini = arma::accu( -h_ % first_sum + second_sum );
        
        // clinical version:
        //be_prop( j ) = R::rnorm( beta_prop_me_(j), beta_prop_sd );
        be_prop( j ) = arma::randn( arma::distr_param( beta_prop_me_(j), beta_prop_sd ) );
        xbeta_prop = xbeta_ - x_.col(j) * be_(j) + x_.col(j) * be_prop(j);
        xbeta_prop.elem( arma::find( xbeta_prop > 700) ).fill( 700. );
        exp_xbeta_prop = arma::exp( xbeta_prop );
        first_sum_prop = arma::sum( matProdVec( ind_r_d_, exp_xbeta_prop ).t(), 1);
        
        h_exp_xbeta_prop_mat = - arma::kron( exp_xbeta_prop, h_.t() );
        h_exp_xbeta_prop_mat.elem( arma::find(h_exp_xbeta_prop_mat > -1.0e-7) ).fill( -1.0e-7 );
        h_exp_xbeta_prop_mat = 1.0 - arma::exp(h_exp_xbeta_prop_mat);
        second_sum_prop = arma::sum( ( arma::log(h_exp_xbeta_prop_mat) % ind_d_ ).t(), 1 );
        
        loglh_prop = arma::accu( - h_ % first_sum_prop + second_sum_prop );
        /*logprior_prop = R::dnorm( be_prop(j), 0.0, sd_be_(j), true);
        logprior_ini = R::dnorm( be_(j), 0.0, sd_be_(j), true);
        logprop_prop = R::dnorm( be_prop(j), beta_prop_me_(j), beta_prop_sd, true);
        logprop_ini = R::dnorm( be_(j), beta_prop_me_(j), beta_prop_sd, true);*/
        logprior_prop = arma::log_normpdf( be_prop(j), 0.0, sd_be_(j) );
        logprior_ini = arma::log_normpdf( be_(j), 0.0, sd_be_(j) );
        logprop_prop = arma::log_normpdf( be_prop(j), beta_prop_me_(j), beta_prop_sd );
        logprop_ini = arma::log_normpdf( be_(j), beta_prop_me_(j), beta_prop_sd );
        logR = loglh_prop - loglh_ini + logprior_prop - logprior_ini + logprop_ini - logprop_prop;
        
        //if( log( R::runif(0., 1.) ) < logR )
        if( log( arma::randu() ) < logR )
        {
            be_(j) = be_prop(j);
            xbeta_ = xbeta_prop;
            sampleRPg_accept_(j)++;
        }
    }
}

// main function
// (i) import data and parameters; (ii) MCMC algorithm; (iii) export estimates
Rcpp::List drive( const std::string& dataFile, const int p, const int q, const std::string& hyperParFile, const std::string& outFilePath,
                 const arma::vec ini_beta, const arma::vec ini_tauSq, const arma::vec ini_h, const arma::uvec groupInd,
          unsigned int nIter, unsigned int nChains, unsigned int thin, bool rw)
{
    
    // set random seed
    #ifdef _OPENMP
        omp_init_lock(&RNGlock);  // init RNG lock for the parallel part
        //arma::arma_rng::set_seed(123);
    #endif
    
    // ###########################################################
    // ## Read Arguments and Data
    // ###########################################################
    
    //Rcout << std::fixed << std::setprecision(11);
    
    Chain_Data chainData; // this initialises the pointers and the strings to ""
    chainData.outFilePath = outFilePath;
    
    // read Data and format into usables
    //Rcout << "Reading input files ... " <<  "\n";
    
    Utils::formatData(dataFile, chainData.survData );
    Utils::readHyperPar(hyperParFile, chainData );
    
    // declare all the data-related variables
    double eta0 = chainData.eta0;
    double kappa0 = chainData.kappa0;
    double c0 = chainData.c0;
    double r = chainData.r;
    double delta = chainData.delta;
    
    arma::vec dataTime = chainData.survData.data->col( 0 );
    arma::vec dataDi = chainData.survData.data->col( 1 );
    arma::vec s0 = {0.}; //event times that are actually observed
    s0(0) = 2. * max( dataTime ) - max( dataTime.elem( arma::find(dataTime != max(dataTime)) ) );
    arma::vec s = arma::join_cols( arma::sort( arma::unique( dataTime.elem( arma::find(dataDi==1.) ) ) ), s0 );
    unsigned int J = s.n_elem;
    arma::uvec groupNo = arma::unique( groupInd );
    unsigned int K = groupNo.n_elem;
    
    ind_r_d = ind_r = ind_d = arma::zeros<arma::mat>( chainData.survData.data->n_rows, J );
    d = arma::sum( ind_d.t(), 1 );
    settingInterval_cpp( dataTime, dataDi, s, J, ind_d, ind_r, ind_r_d, d );
    
    arma::vec beta_prop_me;
    be = beta_prop_me = ini_beta;
    arma::mat beta_p = arma::zeros<arma::mat>( (int)(nIter/thin) + 1, p + q );
    beta_p.row( 0 ) = ini_beta.t();
    
    double lambdaSq;
    lambdaSq = chainData.lambdaSq;
    double sigmaSq;
    sigmaSq = chainData.sigmaSq;
    arma::vec tauSq = ini_tauSq;
    h = ini_h;
    
  // tausq: only for genomic variables
    arma::vec tauSq_exp = arma::zeros<arma::vec>( p );
    for( unsigned int i=0; i<K; ++i )
    {
       // groupNoLocation = arma::find( groupInd == groupNo(i) );
        tauSq_exp.elem( arma::find( groupInd == groupNo(i) ) ).fill( tauSq(i) );
    }
   
    double beta_prop_sd = sqrt(chainData.beta_prop_var);
    double beta_clin_sd = sqrt(chainData.beta_clin_var);
    arma::vec sd_bePart2 = arma::zeros<arma::vec>( q );
    sd_bePart2.fill( beta_clin_sd );
    arma::vec sd_be = arma::join_cols( sqrt(sigmaSq * tauSq_exp), sd_bePart2 );
    
    //arma::uvec xId = arma::linspace<arma::uvec>(2, 1+p+q, p+q)
    arma::mat x = chainData.survData.data->cols( arma::linspace<arma::uvec>(2, 1+p+q, p+q) );
    arma::vec xbeta = x * ini_beta;
    
    arma::vec be_normSq = arma::zeros<arma::vec>( K );
    for( unsigned int i=0; i<K; ++i )
        be_normSq(i) = arma::accu( ini_beta.elem( arma::find(groupInd == groupNo(i)) ) % ini_beta.elem( arma::find(groupInd == groupNo(i)) ) );
    
    if( !any( be_normSq ) )
        be_normSq(0) = 0.1;
    
    arma::vec alpha0 = arma::zeros<arma::vec>( 1 + J );
    //H_star(j) = eta0 * arma::pow( s, kappa0 );
    alpha0.subvec( 1, J ) = c0 * eta0 * arma::pow( s, kappa0 );
    arma::vec hPriorSh = arma::diff( alpha0 );
    
    // save mcmc intermediate estimates
    arma::mat h_p = arma::zeros<arma::mat>( (int)(nIter/thin) + 1, h.n_elem );
    h_p.row( 0 ) = h.t();
    arma::mat tauSq_p = arma::zeros<arma::mat>( (int)(nIter/thin) + 1, tauSq.n_elem );
    tauSq_p.row( 0 ) = tauSq.t();
    
    arma::vec sigmaSq_p = arma::zeros<arma::vec>( (int)(nIter/thin) + 1 );
    arma::vec lambdaSq_p = sigmaSq_p;
    sampleRPc_accept = arma::zeros<arma::uvec>( q );
    sampleRPg_accept = arma::zeros<arma::uvec>( p );
    
    
    // ###########################################################
    // ## MCMC loop
    // ###########################################################
    
    //unsigned int tick = 100; // how many iter for each print?
    unsigned int j = 0; // count thinned results
    const int cTotalLength = 50;
    Rcout << "Running MCMC iterations ...\n";
    
    for( unsigned int M=0; M<nIter; ++M )
    {
        
        if( M % 10 == 0 || M == (nIter - 1) )
            Rcout << "\r[" <<   //'\r' aka carriage return should move printer's cursor back at the beginning of the current line
            std::string(cTotalLength * (M+1.) / nIter, '#') <<        //printing filled part
            std::string(cTotalLength * (1. - (M+1.) / nIter), '-') <<  //printing empty part
            "] " << int( (M+1.) / nIter * 100.0 ) << "%\r"; //printing percentage

        //Updating regression coefficients and hyperparameters
        
        if( q > 0 )
            updateRP_clinical_cpp( p, q, x, ind_r, ind_d, ind_r_d, J, beta_prop_me, beta_prop_sd, xbeta, be, h, sd_be, sampleRPc_accept );
        
        if( rw )
        {
            updateRP_genomic_rw_cpp( p, x, ind_r, ind_d, ind_r_d, J, beta_prop_me, beta_prop_sd, xbeta, be, h, sd_be, sampleRPg_accept );
        } else {
            updateRP_genomic_cpp( p, x, ind_r, ind_d, ind_r_d, J, xbeta, be, h, sd_be, sampleRPg_accept );
        }
        
//        if( q > 0 )
//        {
            for( unsigned int i=0; i<K; ++i )
                be_normSq(i) = arma::accu( be.elem( arma::find(groupInd == groupNo(i)) ) % be.elem( arma::find(groupInd == groupNo(i)) ) );
//        } else {
//            be_normSq = be % be;
//        }
        
        h = updateBH_cpp( ind_r_d, hPriorSh, d, c0, J, xbeta );
            
        tauSq = updateTau_GL_cpp( lambdaSq, sigmaSq, be_normSq );

        sigmaSq = updateSigma_GL_cpp( p, be_normSq, tauSq);
        lambdaSq = updateLambda_GL_cpp( p, K, r, delta, tauSq);
        
//        if( q > 0 )
//        {
            for( unsigned int i=0; i<K; ++i )
                tauSq_exp.elem( arma::find(groupInd == groupNo(i)) ).fill( tauSq(i) );
            
            sd_bePart2.fill( beta_clin_sd );
            sd_be = arma::join_cols( sqrt(sigmaSq * tauSq_exp), sd_bePart2 );
//        } else {
//            tauSq_exp = tauSq;
//            sd_be = sqrt(sigmaSq * tauSq_exp);
//        }

        // Save all results
        if( M % thin == 0 )
        {
            ++j;
            beta_p.row( j )  = be.t();
            h_p.row( j )  = h.t();
            tauSq_p.row( j )  = tauSq.t();
            sigmaSq_p( j ) = sigmaSq;
            lambdaSq_p( j ) = lambdaSq;
        }
        
        // adaptve jumping rule
        if( j > 20 )
            for( unsigned int jj=0; jj<p+q; ++jj )
                if( ini_beta( jj ) == beta_p( j - 20, jj ) )
                    beta_prop_me( jj ) = be( jj );
        
    }
    
    arma::uvec accept_rate = join_cols(sampleRPg_accept, sampleRPc_accept);// / (double)(nIter);
    
    // Exit
    Rcout << "\nDONE, exiting! \n";
    
    return Rcpp::List::create(Rcpp::Named("beta.p") =  beta_p,
                              Rcpp::Named("h.p") = h_p,
                              Rcpp::Named("tauSq.p") = tauSq_p,
                              Rcpp::Named("sigmaSq.p") = sigmaSq_p,
                              Rcpp::Named("lambdaSq.p") = lambdaSq_p,
                              Rcpp::Named("accept.rate") = accept_rate
                              );
    
}

// [[Rcpp::export]]
Rcpp::List psbcSpeedUp_internal( const std::string& dataFile, const int p, const int q, const std::string& hyperParFile, const std::string& outFilePath,
                                const arma::vec ini_beta, const arma::vec ini_tauSq, const arma::vec ini_h, const arma::uvec groupInd, unsigned int nIter, unsigned int nChains, unsigned int thin, bool rw )
{
  //int status {1};
    Rcpp::List beta_mcmc;
    try
    {
        // status = drive(...);
        beta_mcmc =  drive( dataFile, p, q, hyperParFile, outFilePath, ini_beta, ini_tauSq, ini_h, groupInd, nIter, nChains, thin, rw );
    }
    catch( const std::exception& e )
    {
        Rcout << e.what() << '\n'; // we can use Rcerr here because we're reaching here from R for sure
    }
    
    // return status;
    return beta_mcmc;
}
