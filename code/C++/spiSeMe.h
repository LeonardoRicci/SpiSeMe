// SpiSeMe PACKAGE (C++ implementation)
//
//	This header defines the functions provided by the SpiSeMe package.
//	Default values for optional parameters are declared here.

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstring>
#include <vector>
#include <algorithm>
#include <iostream>

enum spiSeMe_return_code : int {
	SSM_SUCCESS = 0,
	SSM_BAD_ARGUMENT = 1,
	SSM_RUNTIME_ERR = 2
};


spiSeMe_return_code spiSeMe_surrogate_jodi(std::vector < std::vector<double> > &, const std::vector<double> &, unsigned int = 1, bool = true);
spiSeMe_return_code spiSeMe_surrogate_iaaft(std::vector < std::vector <double> > &, const std::vector <double> &, std::string = "distribution", unsigned int = 1, bool = true);
spiSeMe_return_code spiSeMe_surrogate_sa(std::vector< std::vector<double> > &, const std::vector<double> &, double, double, double = 0.9, double = 0.1, double = -1.0, std::string = "max", unsigned int = 0, unsigned int = 0, unsigned int = 1, bool = true);
spiSeMe_return_code spiSeMe_surrogate_dither(std::vector < std::vector<double> > &, const std::vector<double> &, std::string = "uniform", double = -1.0, unsigned int = 1, bool = true);

spiSeMe_return_code spiSeMe_event_autocorrelation(std::vector<double> &, std::vector<double> &, const std::vector<double> &, double, double);
spiSeMe_return_code spiSeMe_event_cross_correlation(std::vector<double> &, std::vector<double> &, const std::vector<double> &, const std::vector<double> &, double, double);
spiSeMe_return_code spiSeMe_distribution_test(double &, double &, const std::vector<double> &, const std::vector<double> &);
