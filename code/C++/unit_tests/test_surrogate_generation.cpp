#ifndef SPISEME_HEADER
	#define SPISEME_HEADER
	#include "spiSeMe.h"
#endif
#include <fstream>
#include <sstream>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Test surrogate generation functions"){

	// Load data
	std::ifstream input_file("test_data_A.dat", std::ifstream::in);
	double			temp_iei;
	std::string 		line;
	std::istringstream	linestream;
	std::vector<double>	data_A;
	while (std::getline(input_file, line)) {
		linestream.str(line);
		linestream >> temp_iei;
		data_A.push_back(temp_iei);
		linestream.clear();
	}
	input_file.close();

	std::vector< std::vector<double> > surrogates;
	std::vector<double>	data_Z(1, 0.0);
	spiSeMe_return_code	err;
	unsigned int M = 3;

	SECTION("Surrogate test JODI (arg)") {
		err = spiSeMe_surrogate_jodi(surrogates, data_Z, 1, true);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_jodi(surrogates, data_A, 0, true);
        	REQUIRE(err == SSM_BAD_ARGUMENT);
	}

	SECTION("Surrogate test JODI (run)") {
		err = spiSeMe_surrogate_jodi(surrogates, data_A, M, false);
		REQUIRE(err == SSM_SUCCESS);
		REQUIRE(surrogates.size() == M);
		double p, d;
		err = spiSeMe_distribution_test(p, d, data_A, surrogates[0]);
		REQUIRE(p == 1.0);
	}

	SECTION("Surrogate test IAAFT (arg)") {
		err = spiSeMe_surrogate_iaaft(surrogates, data_Z, "distribution", 1, true);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_iaaft(surrogates, data_A, "distribution", 0, true);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_iaaft(surrogates, data_A, "invalid", M, true);
		REQUIRE(err == SSM_BAD_ARGUMENT);
	}

	SECTION("Surrogate test IAAFT (run)") {
		err = spiSeMe_surrogate_iaaft(surrogates, data_A, "distribution", M, false);
		REQUIRE(err == SSM_SUCCESS);
		REQUIRE(surrogates.size() == M);
		double p, d;
		err = spiSeMe_distribution_test(p, d, data_A, surrogates[0]);
		REQUIRE(p == 1.0);
	}

	SECTION("Surrogate test SA (arg)") {
		err = spiSeMe_surrogate_sa(surrogates, data_Z, 0.1, 1.0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_sa(surrogates, data_A, -0.1, 1.0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_sa(surrogates, data_A, 0.1, -1.0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_sa(surrogates, data_A, 1.0, 0.1);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_sa(surrogates, data_A, 0.1, 1.0, 1.1);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_sa(surrogates, data_A, 0.1, 1.0, -1.1);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_sa(surrogates, data_A, 0.1, 1.0, 0.9, -1.0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_sa(surrogates, data_A, 0.1, 1.0, 0.9, 0.1, 1.0, "invalid");
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_sa(surrogates, data_A, 0.1, 1.0, 0.9, 0.1, 1.0, "max", 100, 50, 0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
	}

	SECTION("Surrogate test SA (run)") {
		err = spiSeMe_surrogate_sa(surrogates, data_A, 0.1, 1.0, 0.9, 0.1, 1.0, "max", 100, 50, M, false);
		REQUIRE(err == SSM_SUCCESS);
		REQUIRE(surrogates.size() == M);
		double p, d;
		err = spiSeMe_distribution_test(p, d, data_A, surrogates[0]);
		REQUIRE(p == 1.0);
	}

	SECTION("Surrogate test Dither (arg)") {
		err = spiSeMe_surrogate_dither(surrogates, data_Z);
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_dither(surrogates, data_A, "invalid");
		REQUIRE(err == SSM_BAD_ARGUMENT);
		err = spiSeMe_surrogate_dither(surrogates, data_A, "uniform", 0.1, 0);
		REQUIRE(err == SSM_BAD_ARGUMENT);
	}

	SECTION("Surrogate test Dither (run)") {
		err = spiSeMe_surrogate_dither(surrogates, data_A, "uniform", 0.1, M, false);
		REQUIRE(err == SSM_SUCCESS);
		REQUIRE(surrogates.size() == M);
	}

}
