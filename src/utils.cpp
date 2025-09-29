#include "utils.h"


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

}
