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

	void readHyperPar(const std::string &hyperParFile, Chain_Data &chainData)
	{
		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_file(hyperParFile.c_str());

		if (result)
		{
			// We'll use Xpaths to get variables
			// the results of empty queries will be a nan after the evaluate_number() and thus we'll preserve NANs for non specified parameters
			pugi::xpath_query query_eta0("/hyperparameters/eta0");
			pugi::xpath_query query_kappa0("/hyperparameters/kappa0");
			pugi::xpath_query query_c0("/hyperparameters/c0");
			pugi::xpath_query query_r("/hyperparameters/r");
			pugi::xpath_query query_delta("/hyperparameters/delta");
			pugi::xpath_query query_beta_prop_var("/hyperparameters/beta_prop_var");
			pugi::xpath_query query_beta_clin_var("/hyperparameters/beta_clin_var");

			pugi::xpath_query query_lambdaSq("/hyperparameters/lambdaSq");
			pugi::xpath_query query_rate("/hyperparameters/rate");
			pugi::xpath_query query_sigmaSq("/hyperparameters/sigmaSq");

			// get results
			chainData.eta0 = query_eta0.evaluate_number(doc);
			chainData.kappa0 = query_kappa0.evaluate_number(doc);
			chainData.c0 = query_c0.evaluate_number(doc);
			chainData.r = query_r.evaluate_number(doc);
			chainData.delta = query_delta.evaluate_number(doc);
			chainData.beta_prop_var = query_beta_prop_var.evaluate_number(doc);
			chainData.beta_clin_var = query_beta_clin_var.evaluate_number(doc);

			chainData.lambdaSq = query_lambdaSq.evaluate_number(doc);
			chainData.rate = query_rate.evaluate_number(doc);
			chainData.sigmaSq = query_sigmaSq.evaluate_number(doc);

			std::vector<std::string> valid_top_level = {"hyperparameters"}; // ,"model","chain"}; ?
																																			// std::vector<std::string> valid_hyperpar = {"eta0","kappa0","c0","r","delta","s","groupInd","beta_prop_var","beta_clin_var", "lambdaSq","rate","beta_ini","sigmaSq","tauSq","h"};
			std::vector<std::string> valid_hyperpar = {"eta0", "kappa0", "c0", "r", "delta", "beta_prop_var", "beta_clin_var", "lambdaSq", "rate", "sigmaSq"};

			for (pugi::xml_node node = doc.first_child(); node; node = node.next_sibling())
			{
				if (std::find(valid_top_level.begin(), valid_top_level.end(), node.name()) == valid_top_level.end())

					Rcout << "\n\n\tWarning: " << node.name() << " not recognised as a valid top level node - only 'hyperparameters' is valid" << '\n'; // ,model and chain
			}

			for (pugi::xml_node node = doc.child("hyperparameters").first_child(); node; node = node.next_sibling())
			{
				if (std::find(valid_hyperpar.begin(), valid_hyperpar.end(), node.name()) == valid_hyperpar.end())
				{
					Rcout << "\n\n\tWARNING: " << node.name() << " was not recognised as a valid hyperaparameter" << '\n';
					Rcout << "\t" << node.name() << ": " << node.child_value() << " disregarded .. " << '\n';
					Rcout << "\tValid hyperparameters are: \n\t";
					for (auto &v : valid_hyperpar)
						Rcout << v << ", ";
					Rcout << " --- see the documentation for more details " << '\n'
								<< '\n';
				}
			}
		}
		else
		{
			Rcout << '\n'
						<< "No hyperparameter input file was given (or wrong format detected), so default values will be used." << '\n';
		}
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
