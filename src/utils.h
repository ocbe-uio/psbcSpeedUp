#ifndef UTILS
#define UTILS

#ifdef CCODE
    #include <iostream>
    #include <armadillo>
#else
    #include <RcppArmadillo.h>
#endif

#include <memory>
#include <string>
#include <cmath>
#include <limits>

#include "pugixml.h"

namespace Utils{

	struct Surv_Data
	{
		std::shared_ptr<arma::mat> data;
		//unsigned int nObservations, nPredictors, nVSPredictors, nFixedPredictors;
        
		//std::shared_ptr<arma::umat> missingDataArrayIdx;
		//std::shared_ptr<arma::uvec> completeCases;

		Surv_Data() // use this constructor to instanciate all the object at creation (to be sure pointers point to *something*)
		{
			data = std::make_shared<arma::mat>();

			//missingDataArrayIdx = std::make_shared<arma::umat>();
			//completeCases = std::make_shared<arma::uvec>();
		}
	};

	struct Chain_Data
	{
		// Data
		Utils::Surv_Data survData;

		// Misc MCMC quantities
		unsigned int nChains = 1, nIter = 10, burnin = 0;
		// hyperparameters

        double eta0		= std::nan("0");
		double kappa0 	= std::nan("0");
		double c0       = std::nan("0");
		double r		= std::nan("0");
		double delta 	= std::nan("0");
        /*double eta0 = 1.;
        double kappa0 = 1.;
        double c0 = 1.;
        double r = 1.;
        double delta = 1.;*/
        
        // init for some variables
        double lambdaSq = 1.;
        double rate = lambdaSq / 2.;
        double sigmaSq = 1.;
        //arma::vec beta_ini;
        //arma::vec tauSq;
        //arma::vec h;
        
        // other parameters
        double beta_prop_var = 1.;
        double beta_clin_var = 1.;
        
		// file names and paths
		std::string filePrefix , outFilePath;

		// outputs
		//bool output_gamma, output_beta, output_sigmaRho,
		//	output_Gy, output_pi, output_tail, output_model_size, output_CPO, output_model_visit;
        
	};

//    struct Initial_Data
//    {
//        // Data
//        Utils::Ini_Data iniData;
//
//        // init for some variables
//        double lambdaSq = 1;
//        double rate = lambdaSq / 2.;
//        arma::vec beta_ini;
//        arma::vec sigmaSq;
//        arma::vec tauSq;
//        arma::vec h;
//
//    };

	class badFile : public std::exception
	{
		const char * what () const throw ()
		{
			return "The file is either missing or in a wrong format, make sure you're feeding plaintext files.";
		}
	};

	class badRead : public std::exception
	{
		const char * what () const throw ()
		{
			return "Unknown error: Something went wrong while reading the files. Check your input.";
		}
	};
	

	bool readData(const std::string& dataFileName, std::shared_ptr<arma::mat> data);
    
	/* Computes the set-difference from two vectors of indexes */
	arma::uvec arma_setdiff_idx(const arma::uvec& x, const arma::uvec& y);

	//void initMissingData(std::shared_ptr<arma::mat> data, std::shared_ptr<arma::umat> missingDataArrayIndexes, std::shared_ptr<arma::uvec> completeCases, bool print=false );

	void formatData(const std::string& dataFileName, Surv_Data& survData );

	void readHyperPar(const std::string& hyperParFile, Chain_Data& chainData );
    //void readInitialPar(const std::string& initialParFile, Initial_Data& initialData );

	template <typename T> int sgn(T val)
	{
		return (T(0) < val) - (val < T(0));
	}

	inline double internalRound( double val )
	{
		if( val < 0 ) return ceil(val - 0.5);
		return floor(val + 0.5);
	}

	inline double round( double val , unsigned int decimals )
	{
		return internalRound( val * std::pow(10.,decimals) ) / std::pow(10.,decimals) ;
	}

	double logspace_add(const arma::vec& logv);
	double logspace_add(double a,double b);

	arma::uvec nonZeroLocations_row( arma::sp_umat X);  // if you pass a row subview
	arma::uvec nonZeroLocations_col( arma::sp_umat X); // if you pass a col subview

	// arma::uvec nonZeroLocations_col( arma::sp_umat X, unsigned int column);
	// arma::uvec nonZeroLocations_row( arma::sp_umat X, unsigned int row);
		
}

#endif
