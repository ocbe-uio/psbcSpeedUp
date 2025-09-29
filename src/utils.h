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


namespace Utils
{

	struct Surv_Data
	{
		std::shared_ptr<arma::mat> data;
		// unsigned int nObservations, nPredictors, nVSPredictors, nFixedPredictors;

		// std::shared_ptr<arma::umat> missingDataArrayIdx;
		// std::shared_ptr<arma::uvec> completeCases;

		Surv_Data() // use this constructor to instanciate all the object at creation (to be sure pointers point to *something*)
		{
			data = std::make_shared<arma::mat>();

			// missingDataArrayIdx = std::make_shared<arma::umat>();
			// completeCases = std::make_shared<arma::uvec>();
		}
	};

	struct Chain_Data
	{
		// Data
		Utils::Surv_Data survData;

	};

	class badFile : public std::exception
	{
		const char *what() const throw()
		{
			return "The file is either missing or in a wrong format, make sure you're feeding plaintext files.";
		}
	};

	bool readData(const std::string &dataFileName, std::shared_ptr<arma::mat> data);

	void formatData(const std::string &dataFileName, Surv_Data &survData);

	template <typename T>
	int sgn(T val)
	{
		return (T(0) < val) - (val < T(0));
	}

}

#endif
