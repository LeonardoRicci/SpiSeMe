#ifndef SPISEME_HEADER
	#define SPISEME_HEADER
	#include "spiSeMe.h"
#endif
#include <fstream>
#include <sstream>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Test auxiliary functions"){

	// Load data
	std::ifstream input_file_A("test_data_A.dat", std::ifstream::in);
	std::ifstream input_file_B("test_data_B.dat", std::ifstream::in);
	double			temp_iei;
	std::string 		line;
	std::istringstream	linestream;
	std::vector<double>	data_A;
	std::vector<double>	data_B;
	while (std::getline(input_file_A, line)) {
		linestream.str(line);
		linestream >> temp_iei;
		data_A.push_back(temp_iei);
		linestream.clear();
	}
	input_file_A.close();
	while (std::getline(input_file_B, line)) {
		linestream.str(line);
		linestream >> temp_iei;
		data_B.push_back(temp_iei);
		linestream.clear();
	}
	input_file_B.close();

	std::vector<double>	A;
	std::vector<double>	C;
	std::vector<double>	lags;
	std::vector<double>	data_Z(1, 0.0);
	spiSeMe_return_code	err;

	SECTION("Auxiliary test KS-test (arg)") {
		double p, d;
		err = spiSeMe_distribution_test(p, d, data_A, data_Z);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_distribution_test(p, d, data_Z, data_A);
		REQUIRE(err == SSM_BAD_ARGUMENT);
	}

	SECTION("Auxiliary test KS-test (run)") {
		double p, d;
		err = spiSeMe_distribution_test(p, d, data_A, data_B);
		REQUIRE(err == SSM_SUCCESS);
		REQUIRE(fabs(p-0.5) <= 0.1);
		REQUIRE(fabs(d-0.04) <= 0.01);
	}

	SECTION("Auxiliary test Autocorrelation (arg)") {
		err = spiSeMe_event_autocorrelation(A, lags, data_Z, 0.1, 1.0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_event_autocorrelation(A, lags, data_A, -0.1, 1.0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_event_autocorrelation(A, lags, data_A, 0.1, -1.0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_event_autocorrelation(A, lags, data_A, 1.0, 0.1);
		REQUIRE(err == SSM_BAD_ARGUMENT);
	}

	SECTION("Auxiliary test Autocorrelation (run)") {
		err = spiSeMe_event_autocorrelation(A, lags, data_A, 0.1, 1.0);
		REQUIRE(err == SSM_SUCCESS);
		REQUIRE(A.size() == 10);
		REQUIRE(A.size() == lags.size());
	}

	SECTION("Auxiliary test Cross-Correlation (arg)") {
		err = spiSeMe_event_cross_correlation(C, lags, data_Z, data_A, 0.1, 1.0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_event_cross_correlation(C, lags, data_A, data_Z, 0.1, 1.0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_event_cross_correlation(C, lags, data_A, data_B, -0.1, 1.0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_event_cross_correlation(C, lags, data_A, data_B, 0.1, -1.0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_event_cross_correlation(C, lags, data_A, data_B, 1.0, 0.1);
		REQUIRE(err == SSM_BAD_ARGUMENT);
	}

	SECTION("Auxiliary test Cross-Correlation (run)") {
		err = spiSeMe_event_cross_correlation(C, lags, data_A, data_B, 0.1, 1.0);
		REQUIRE(err == SSM_SUCCESS);
		REQUIRE(C.size() == 21);
		REQUIRE(C.size() == lags.size());
	}

}
