// SPISEME_SURROGATE_IAAFT
//
//	Generates, by means of the Iterative Amplitude Adjusted Fourier Transform
//	(IAAFT) algorithm, surrogate IEI sequences corresponding to the original
//	sequence iei_sequence.
//
//	'exactly_preserve' Specifies which function has to be preserved
//			exactly by the surrogate sequence. Can either be
//			"distribution" (default) or "spectrum".
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
//	The IAAFT (Iterative Amplitude-Adjusted Fourier Transform) algorithm
//	for surrogate generation was originally proposed by T. Schreiber
//	and A. Schmitz in Phys. Rev. Lett. 77 (1996), 635,
//	doi:10.1103/PhysRevLett.77.635
//
//	A comparative study of surrogate generation algorithms can be found in
//	Chaos 29 (2019), 121102, doi:10.1063/1.5138250
//
//	Citing this source is highly appreciated. Thank you.

#ifndef SPISEME_HEADER
	#define SPISEME_HEADER
	#include "spiSeMe.h"
#endif
#include "gsl/gsl_fft_real.h"
#include "gsl/gsl_fft_halfcomplex.h"
#include <gsl/gsl_rng.h>

struct pair_index_value {
	unsigned int	index;
	double		value;
};

bool sort_pairs_by_value(pair_index_value a, pair_index_value b) { return (a.value < b.value); }
bool sort_pairs_by_index(pair_index_value a, pair_index_value b) { return (a.index < b.index); }
void assess_original_spectrum(std::vector<double> &, const std::vector<double> &);
void adjust_spectrum(std::vector<double> &, const std::vector<double> &, const std::vector<double> &);
void adjust_distribution(std::vector<double> &, const std::vector<double> &);
double check_sequence_convergence(const std::vector<double> &, const std::vector<double> &);

spiSeMe_return_code spiSeMe_surrogate_iaaft(std::vector < std::vector <double> > & iei_surrogates, const std::vector <double> & iei_sequence, std::string exactly_preserve, unsigned int M, bool verbose)
{
	// --- Input parsing & validation
	unsigned int N = iei_sequence.size();
	if (N < 2) {
		std::cerr << "ERROR (in spiSeMe_surrogate_iaaft): input iei_sequence is not an actual vector (length < 2).\n";
		return SSM_BAD_ARGUMENT;
	}
	if ((strcmp(exactly_preserve.c_str(), "distribution")) && (strcmp(exactly_preserve.c_str(), "spectrum"))) {
		std::cerr << "ERROR (in spiSeMe_surrogate_iaaft): \n";
		return SSM_BAD_ARGUMENT;
	}
	if (M < 1) {
		std::cerr << "ERROR (in spiSeMe_surrogate_iaaft): function argument M must be a positive integer.\n";
		return SSM_BAD_ARGUMENT;
	}

	// --- Provide some information
	if (verbose) {
		std::cerr << "\n### Starting IAAFT routine ###\n";
		std::cerr << "###\tExactly preserving: " << exactly_preserve << "\n";
	}

	// --- Assess spectrum & distribution of original sequence
	std::vector <double>	spectrum_magnitude_original;
	std::vector <double>	distribution_original(iei_sequence);
	std::sort(distribution_original.begin(), distribution_original.end());
	assess_original_spectrum(spectrum_magnitude_original, iei_sequence);

	// --- Initialize output array (each of the M generated surrogates is a vector)
	iei_surrogates.clear();
	std::vector<double> sequence_elements;
	std::vector<double> temp_iei_surrogate;
	std::vector<double> iei_surrogate_spectrum_adjusted;

	// --- Initialize random engine
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, time(NULL));

	//--- Generate surrogates: iterate M times
       for (int iter = 0; iter < M; iter++) {
	       if (verbose)
		       std::cerr << "# Surrogate number " << iter+1 << " out of " << M << ".\n";

	       // --- Initialize surrogate ranks generation
	       bool run_iterations = true;
	       int idx;
	       double discrepancy;
	       temp_iei_surrogate.clear();
	       for (int i = 0; i < N; i++) {
		       sequence_elements.push_back(iei_sequence[i]);
	       }
	       for (int i = 0; i < N; i++) {
		       idx = gsl_rng_uniform_int(r, sequence_elements.size());
		       temp_iei_surrogate.push_back(sequence_elements[idx]);
		       sequence_elements.erase(sequence_elements.begin() + idx);
	       }

		// --- Iterative algorithm (switched depending on target function to be matched exactly.
		int n_iter = 0;
		while (run_iterations) {
			n_iter = n_iter + 1;

			adjust_spectrum(iei_surrogate_spectrum_adjusted, temp_iei_surrogate, spectrum_magnitude_original);
			adjust_distribution(iei_surrogate_spectrum_adjusted, distribution_original);

			discrepancy = check_sequence_convergence(iei_surrogate_spectrum_adjusted, temp_iei_surrogate);
			if (discrepancy == 0) {
				if (!strcmp(exactly_preserve.c_str(), "spectrum")) {
					adjust_spectrum(temp_iei_surrogate, iei_surrogate_spectrum_adjusted, spectrum_magnitude_original);
				}
				break;
			} else if (discrepancy < 0) {
				std::cerr << "ERROR (in spiSeMe_surrogate_iaaft): error during convergence check.\n";
				return SSM_RUNTIME_ERR;
			}
				std::copy(iei_surrogate_spectrum_adjusted.begin(), iei_surrogate_spectrum_adjusted.end(), temp_iei_surrogate.begin());
		}

		if (verbose)
			std::cerr << "Process converged after " << n_iter << " iterations.\n\n";

		// --- Append sequence to the set of those already generated
		iei_surrogates.push_back(temp_iei_surrogate);
       }

	gsl_rng_free(r);

	return SSM_SUCCESS;
}

void assess_original_spectrum(std::vector<double> & spectrum_magnitude_original, const std::vector<double> & sequence)
{
	unsigned int N = sequence.size();
	double * data = (double *) malloc(sizeof(double) * N);
	for (int i = 0; i < N; i++) {
		data[i] = sequence[i];
	}
	gsl_fft_real_wavetable * wt = gsl_fft_real_wavetable_alloc(N);
	gsl_fft_real_workspace * workspace = gsl_fft_real_workspace_alloc(N);

	gsl_fft_real_transform(data, 1, N, wt, workspace);
	spectrum_magnitude_original.clear();
	spectrum_magnitude_original.push_back(data[0]);		// Zero-freq., purely real
	for (int i = 1; i < N; i++) {
		if ((i < N-1) || (N%2 != 0)) {			// Almost all cases (all cases if N odd)
			spectrum_magnitude_original.push_back(sqrt(data[i]*data[i] + data[i+1]*data[i+1]));
			i++;
		} else {					// If N even, last coeff is purely real as well
			spectrum_magnitude_original.push_back(data[i]);
		}
	}

	gsl_fft_real_wavetable_free(wt);
	gsl_fft_real_workspace_free(workspace);
	free(data);

	return;
}

void adjust_spectrum(std::vector <double> & adjusted_sequence, const std::vector<double> & sequence, const std::vector<double> & original_amplitudes)
{
	unsigned int N = sequence.size();
	unsigned int K = original_amplitudes.size();

	double * data = (double *) malloc(sizeof(double) * N);
	for (int i = 0; i < N; i++) {
		data[i] = sequence[i];
	}

	gsl_fft_real_wavetable * wtr = gsl_fft_real_wavetable_alloc(N);
	gsl_fft_real_workspace * workspace = gsl_fft_real_workspace_alloc(N);
	gsl_fft_real_transform(data, 1, N, wtr, workspace);

	double angle;
	for (int i = 0; i < K; i++) {
		if (i == 0) {					// Zero-freq., purely real
			data[0] = original_amplitudes[0];
		} else if ((N%2 != 0) || (i < K - 1)) {		// Almost all cases (all cases if N odd)
			angle = atan2(data[2*i], data[2*(i-1)+1]);
			data[2*(i-1)+1] = original_amplitudes[i] * cos(angle);
			data[2*i] = original_amplitudes[i] * sin(angle);
		} else {
			data[2*(i-1)+1] = original_amplitudes[i];	// If N even, last coeff is purely real as well
		}
	}

	gsl_fft_halfcomplex_wavetable * wthc = gsl_fft_halfcomplex_wavetable_alloc(N);
	gsl_fft_halfcomplex_inverse(data, 1, N, wthc, workspace);
	gsl_fft_real_wavetable_free(wtr);
	gsl_fft_real_workspace_free(workspace);
	gsl_fft_halfcomplex_wavetable_free(wthc);

	adjusted_sequence.clear();
	for (int i = 0; i < N; i++)
		adjusted_sequence.push_back(data[i]);

	free(data);

	return;
}

void adjust_distribution(std::vector<double> & sequence, const std::vector<double> & distribution_original)
{
	unsigned int N = sequence.size();

	std::vector <pair_index_value>	indexed_sequence;
	pair_index_value	temp_pair;
	for (int i = 0; i < N; i++) {
		temp_pair.index = i;
		temp_pair.value = sequence[i];
		indexed_sequence.push_back(temp_pair);
	}

	std::sort(indexed_sequence.begin(), indexed_sequence.end(), sort_pairs_by_value);
	for (int i = 0; i < N; i++) {
		indexed_sequence[i].value = distribution_original[i];
	}
	std::sort(indexed_sequence.begin(), indexed_sequence.end(), sort_pairs_by_index);

	for (int i = 0; i < N; i++) {
		sequence[i] = indexed_sequence[i].value;
	}

	return;
}

double check_sequence_convergence(const std::vector<double> & sequence_a, const std::vector<double> & sequence_b)
{
	// --- Define check of convergence (sequence did not change with respect to last iteration)
	unsigned int N = sequence_a.size();
	if (N != sequence_b.size()) {
		std::cerr << "ERROR (in spiSeMe_surrogate_iaaft): internal error due to size mismatch when checking sequence convergence.\n";
		return -1;
	}

	double difference = 0.0;
	for (int i = 0; i < N; i++) {
		difference += fabs(sequence_a[i] - sequence_b[i]);
	}

	return difference;
}
