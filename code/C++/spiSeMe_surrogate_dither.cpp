// SPISEME_SURROGATE_DITHER
//
//	Generates, by means of dithering, surrogate IEI sequences corresponding
//	to the original sequence iei_sequence.
//
//	Dithering consists in adding random time intervals to the time
//	coordinate of each event in the sequence. Random intervals follow a
//	zero-mean distribution whose shape and width can be set by the
//	user.
//
//	'dither_distribution' Specifies from which distribution have the dithering
//			elements to be drawn. Can either be "uniform"
//			(default), "normal", or "triangular".
//
//	'distribution_parameter' In the case of "uniform" or "triangular", sets
//			the half-width of the distribution from	which dithering
//			elements have to be drawn (see 'dither_distribution'). 
//			In the case of "normal", sets the distribution standard
//			deviation. By default, D is estimated as half the minimum
//			IEI within the original sequence.
//
//	'M'		Number of surrogate sequences to generate.
//			Default is 1. Each of the M surrogate sequences
//			is an element of the output vector<vector<>>, which
//			therefore has outer size M and all inner sizes N,
//			where N is the original sequence length.
//
//	'verbose'	Sets the verbosity of the function. If true (default),
//			all messages are printed on the command line. If
//			false, only critical errors are displayed.
//
//	This function is part of the SpiSeMe package.
//
//
//	REFERENCE:
//	A comparative study of surrogate generation algorithms can be found in
//	Chaos 29 (2019), 121102, doi:10.1063/1.5138250
//
//	Citing this source is highly appreciated. Thank you.

#ifndef SPISEME_HEADER
	#define SPISEME_HEADER
	#include "spiSeMe.h"
#endif
#include <gsl/gsl_randist.h>

spiSeMe_return_code spiSeMe_surrogate_dither(std::vector < std::vector<double> > & iei_surrogates, const std::vector<double> & iei_sequence, std::string dither_distribution, double distribution_parameter, unsigned int M, bool verbose)
{
	// --- Input parsing & validation
	unsigned int N = iei_sequence.size();
	if (N < 2) {
		std::cerr << "ERROR (in spiSeMe_surrogate_dither): input iei_sequence is not an actual vector (length < 2).\n";
		return SSM_BAD_ARGUMENT;
	}
	if (M < 1) {
		std::cerr << "ERROR (in spiSeMe_surrogate_dither): function argument M must be a positive integer.\n";
		return SSM_BAD_ARGUMENT;
	}
	if (strcmp(dither_distribution.c_str(), "uniform") && strcmp(dither_distribution.c_str(), "normal") && strcmp(dither_distribution.c_str(), "triangular")) {
		std::cerr << "ERROR (in spiSeMe_surrogate_dither): function argument <dither_distribution> does not match any of the allowed values (''uniform'', ''normal'', ''triangular''.\n";
		return SSM_BAD_ARGUMENT;
	}
	if (distribution_parameter <= 0) {
		double min_iei = std::numeric_limits<double>::max();
		for (int i = 0; i < N; i++) {
			if (iei_sequence[i] < min_iei)
				min_iei = iei_sequence[i];
		}
		distribution_parameter = min_iei / 2.0;
	}

	// --- Provide some information
	if (verbose)
		std::cerr << "\n### Starting dithering routine ###\n";
		std::cerr << "###\tDither distribution: " << dither_distribution << "\n";
		if (!strcmp(dither_distribution.c_str(), "uniform"))
			std::cerr << "###\t\tU(-D,D), with D = " << distribution_parameter << "\n";
		else if (!strcmp(dither_distribution.c_str(), "normal"))
			std::cerr << "###\t\tN(0,s^2), with s = " << distribution_parameter << "\n";
		else if (!strcmp(dither_distribution.c_str(), "triangular"))
			std::cerr << "###\t\tT(-D,D), with D = " << distribution_parameter << "\n";

	// --- Build arrival times to be dithered
	std::vector<double> original_arrival_times(N+1, 0.0);
	for (int i = 0; i < N; i++)
		original_arrival_times[i+1] = original_arrival_times[i] + iei_sequence[i];

	// --- Initialize output array (each of the M generated surrogates is a column)
	iei_surrogates.clear();
	std::vector<double> temp_surrogate;

	// --- Initialize random engine
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, time(NULL));

	// --- Generate surrogates: iterate M times
	for (int iter = 0; iter < M; iter++) {
		if(verbose)
			std::cerr << "# Surrogate number " << iter+1 << " out of " << M << "\n";

		// --- Initialize surrogate ranks generation
		std::vector<double> surrogate_arrival_times(original_arrival_times);

		// --- Apply dithering to arrival times
		if (!strcmp(dither_distribution.c_str(), "uniform")) {
			for (int i = 0; i < surrogate_arrival_times.size(); i++)
				surrogate_arrival_times[i] += distribution_parameter*(-1.0 + 2.0*gsl_rng_uniform(r));
		} else if (!strcmp(dither_distribution.c_str(), "normal")) {
			for (int i = 0; i < surrogate_arrival_times.size(); i++)
				surrogate_arrival_times[i] += gsl_ran_gaussian(r, distribution_parameter);
		} else if (!strcmp(dither_distribution.c_str(), "triangular")) {
			for (int i = 0; i < surrogate_arrival_times.size(); i++)
				surrogate_arrival_times[i] += distribution_parameter * (gsl_rng_uniform(r) - gsl_rng_uniform(r));
		}

		// --- Transform arrival times back into IEIs
		sort(surrogate_arrival_times.begin(), surrogate_arrival_times.end());
		temp_surrogate.clear();
		for (int i = 0; i < N; i++) {
			temp_surrogate.push_back(surrogate_arrival_times[i+1] - surrogate_arrival_times[i]);
		}

		if (verbose)
			std::cerr << "Process ended.\n\n";

		// --- Append sequence to the set of those already generated
		iei_surrogates.push_back(temp_surrogate);
	}

	gsl_rng_free(r);

	return SSM_SUCCESS;
}
