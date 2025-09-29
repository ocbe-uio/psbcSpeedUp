#include "utils.h"

#ifndef CCODE
using Rcpp::Rcerr;
using Rcpp::Rcout;
#else
#define Rcout std::cout
#define Rcerr std::cerr
#endif

namespace Utils
{

	bool readData(const std::string &dataFileName, std::shared_ptr<arma::mat> data)
	{
		bool status = data->load(dataFileName, arma::raw_ascii);
		if (!status)
			throw badFile();

		// nObservations = data->n_rows;

		return status;
	}
	void formatData(const std::string &dataFileName, Surv_Data &surData)
	{

		readData(dataFileName, surData.data);
	}


	// sgn is defined in the header in order for it to be visible

	double logspace_add(const arma::vec &logv)
	{

		if (logv.is_empty())
			return std::numeric_limits<double>::lowest();
		double m;
		if (logv.has_inf() || logv.has_nan()) // || logv.is_empty()
		{
			return logspace_add(logv.elem(arma::find_finite(logv)));
		}
		else
		{
			m = arma::max(logv);
			return m + std::log((double)arma::as_scalar(arma::sum(arma::exp(-(m - logv)))));
		}
	}

	double logspace_add(double a, double b)
	{

		if (a <= std::numeric_limits<float>::lowest())
			return b;
		if (b <= std::numeric_limits<float>::lowest())
			return a;
		return std::max(a, b) + std::log((double)(1. + std::exp((double)-std::abs((double)(a - b)))));
	}

	arma::uvec nonZeroLocations_col(arma::sp_umat X)
	{
		std::vector<arma::uword> locations;

		for (arma::sp_umat::const_iterator it = X.begin(); it != X.end(); ++it)
		{
			locations.push_back(it.row());
		}

		return arma::uvec(locations);
	}

	arma::uvec nonZeroLocations_row(arma::sp_umat X)
	{
		std::vector<arma::uword> locations;

		for (arma::sp_umat::const_iterator it = X.begin(); it != X.end(); ++it)
		{
			locations.push_back(it.col());
		}

		return arma::uvec(locations);
	}

}
