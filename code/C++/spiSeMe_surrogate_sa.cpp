// SPISEME_SURROGATE_SA
//
//	Generates, by means of the Simulated Annealing (SA) algorithm,
//	surrogate IEI sequences corresponding to the original
//	sequence iei_sequence.
//
//	The evaluation of the cost function relies on the assessment of
//	event autocorrelation. This assessment is implemented by the
//	function spiSeMe_event_autocorrelation within this package. The mandatory
//	parameters 'autocorr_bin_width', 'autocorr_max_lag' correspond to
//	the 'bin_width', 'max_lag' parameters of the spiSeMe_event_autocorrelation
//	function: please refer to its documentation for further details.
//
//	'a'		Cooling factor (in references, alpha). Default 0.9.
//
//	'T'		Starting temperature. Default 0.1.
//
//	'C'		Target cost below which the routine is stopped. If
//			not specified, the function attempts to estimate a
//			reasonable target cost as the standard deviation of the
//			original sequence's autocorrelation. This estimate is not
//			necessarily a good choice for every sequence. Similarly,
//			different metrics (see 'cost_function') will lead to a
//			different accuracy under the default target cost.
//
//	'cost_function'	Sets the metric used to evaluate the cost function:
//			'max': Cost equals the maximim absolute difference
//			between original and surrogate autocorrelations
//			(this is the default setting).
//			'L1': Cost equals the sum over all bins of the absolute
//			differences between original and surrogate
//			autocorrelations.
//			'L2': Cost equals the sum over all bins of the
//			squared differences between original and surrogate
//			autocorrelations.
//
//	'n_total'	Number of swaps to be carried out before cooling.
//			Default is N, where N is the sequence length.
//
//	'n_successful'	Number of accepted swaps to be carried out before
//			cooling. Default is N/2, where N is the sequence length.
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
//	The simulated annealing algorithm for surrogate generation was originally
//	proposed by T. Schreiber in Phys. Rev. Lett. 80 (1998), 2105,
//	doi:10.1103/PhysRevLett.80.2105
//
//	Details on the autocorrelation assessment, as well as a comparative study
//	of surrogate generation algorithms, can be found in
//	Chaos 29 (2019), 121102, doi:10.1063/1.5138250
//
//	Citing this source is highly appreciated. Thank you.
//
//	See also EVENT_AUTOCORRELATION

#ifndef SPISEME_HEADER
	#define SPISEME_HEADER
	#include "spiSeMe.h"
#endif
#include <gsl/gsl_rng.h>

double compute_cost_function(const std::vector<double> &, const std::vector<double> &, unsigned int);

spiSeMe_return_code spiSeMe_surrogate_sa(std::vector< std::vector<double> > & iei_surrogates, const std::vector<double> & iei_sequence, double autocorr_bin_width, double autocorr_max_lag, double a, double T, double C, std::string cost_function, unsigned int n_total, unsigned int n_successful, unsigned int M, bool verbose)
{
	// --- Input parsing & validation
	unsigned int N = iei_sequence.size();
	if (N < 2) {
		std::cerr << "ERROR (in spiSeMe_surrogate_sa): input iei_sequence is not an actual vector (length < 2).\n";
		return SSM_BAD_ARGUMENT;
	}
	if (M < 1) {
		std::cerr << "ERROR (in spiSeMe_surrogate_sa): function argument M must be a positive integer.\n";
		return SSM_BAD_ARGUMENT;
	}

	if (autocorr_bin_width <= 0) {
	 	std::cerr << "ERROR (in spiSeMe_event_autocorrelation): argument 'bin_width' must be positive.\n";
		return SSM_BAD_ARGUMENT;
	}
	if (autocorr_max_lag <= 0) {
	 	std::cerr << "ERROR (in spiSeMe_event_autocorrelation): argument 'max_lag' must be positive.\n";
		return SSM_BAD_ARGUMENT;
	}
	if (autocorr_max_lag <= autocorr_bin_width) {
	 	std::cerr << "ERROR (in spiSeMe_event_autocorrelation): invalid parameters. Bin width must be smaller than maximum lag.\n";
		return SSM_BAD_ARGUMENT;
	}

	unsigned int cost_code;
	if (!strcmp(cost_function.c_str(), "max"))
		cost_code = 0;
	else if (!strcmp(cost_function.c_str(), "L1"))
		cost_code = 1;
	else if (!strcmp(cost_function.c_str(), "L2"))
		cost_code = 2;
	else {
		std::cerr << "ERROR (in spiSeMe_event_autocorrelation): invalid parameter. Function argument <cost_function> does not match any of the allowed values (''max'', ''L1'', ''L2'').\n";
		return SSM_BAD_ARGUMENT;
	}

	double cooling_factor;
	if ((a >= 1.0) || (a <= 0)) {
		std::cerr << "ERROR (in spiSeMe_surrogate_sa): invalid parameter. Cooling factor must be positive and less than unity.\n";
		return SSM_BAD_ARGUMENT;
	} else
		cooling_factor = a;
	double starting_temperature;
	if (T < 0) {
		std::cerr << "ERROR (in spiSeMe_surrogate_sa): invalid parameter. Starting temperature must be non-negative.\n";
		return SSM_BAD_ARGUMENT;
	} else
		starting_temperature = T;
	double target_cost = C;


	// --- Assess defaults for omitted parameters
	if (n_total <= 0)
		n_total = N;
	if (n_successful <= 0)
		n_successful = ceil(N / 2);

	// --- Assess autocorrelation of original sequence
	spiSeMe_return_code err_ac;
	std::vector<double> autocorrelation_original;
	std::vector<double> lags;
	err_ac = spiSeMe_event_autocorrelation(autocorrelation_original, lags, iei_sequence, autocorr_bin_width, autocorr_max_lag);
	if (err_ac != SSM_SUCCESS) {
		std::cerr << "ERROR (in spiSeMe_surrogate_sa): error in autocorrelation evaluation.\n";
		return SSM_RUNTIME_ERR;
	}

	// --- If not assigned, estimate the target cost as the standard deviation of the Autocorrelation
	if (target_cost <= 0) {
		double s = 0.0, m = 0.0;
		for (int i = 0; i < N; i++) {
			s += autocorrelation_original[i] * autocorrelation_original[i];
			m += autocorrelation_original[i];
		}
		target_cost = sqrt((s - m*m/N)/(N-1));
		std::cerr << "WARNING (in spiSeMe_surrogate_sa): automatic setting of target cost might be too small for the routine to converge in a reasonable time.\n";
	}

	// --- Provide some information
	if (verbose) {
		std::cerr << "\n### Starting simulated annealing routine ###\n";
		std::cerr << "###\tTarget cost: " << target_cost << "\n";
		std::cerr << "###\tCooling factor: " << cooling_factor << "\n";
		std::cerr << "###\tStarting T: " << starting_temperature << "\n";
		std::cerr << "###\tSuccessful swaps before cooling: " << n_successful << "\n";
		std::cerr << "###\tTotal swaps before cooling: " << n_total << "\n\n";
	}

	// --- Initialize output array (each of the M generated surrogates is a vector)
	iei_surrogates.clear();
	std::vector<double> surrogate_sequence;
	std::vector<double> sequence_elements;

	// --- Initialize random engine
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, time(NULL));

	// --- Generate surrogates: iterate M times
	for (int iter = 0; iter < M; iter++) {
		if (verbose)
			std::cerr << "# Surrogate number " << iter+1 <<" out of " << M <<".\n\n";

		// --- Starting conditions
		double T = starting_temperature;
		int n_swaps = 0;
		bool run_cooling = true;
		bool run_randomization = true;
		int idx, err_cnvg;
		surrogate_sequence.clear();
		sequence_elements.clear();
		for (int i = 0; i < N; i++) {
			sequence_elements.push_back(iei_sequence[i]);
		}
		for (int i = 0; i < N; i++) {
			idx = gsl_rng_uniform_int(r, sequence_elements.size());
			surrogate_sequence.push_back(sequence_elements[idx]);
			sequence_elements.erase(sequence_elements.begin() + idx);
		}
		std::vector<double> autocorrelation_surrogate;
		err_ac = spiSeMe_event_autocorrelation(autocorrelation_surrogate, lags, surrogate_sequence, autocorr_bin_width, autocorr_max_lag);
		if (err_ac != SSM_SUCCESS) {
			std::cerr << "ERROR (in spiSeMe_surrogate_sa): error in autocorrelation evaluation.\n";
			return SSM_RUNTIME_ERR;
		}
		double previous_cost = compute_cost_function(autocorrelation_surrogate, autocorrelation_original, cost_code);

		// --- Outermost iteration: cooling stages
		int Ns, Nt, idx_a, idx_b;
		double permutation_cost, delta_E, swap_temp;
		while (run_cooling) {

			Ns = 0;
			Nt = 0;

			// --- Innermost iteration: swaps & Metropolis step
			while (run_randomization) {

				idx_a = gsl_rng_uniform_int(r, N);
				idx_b = gsl_rng_uniform_int(r, N);
				swap_temp = surrogate_sequence[idx_b];
				surrogate_sequence[idx_b] = surrogate_sequence[idx_a];
				surrogate_sequence[idx_a] = swap_temp;

				err_ac = spiSeMe_event_autocorrelation(autocorrelation_surrogate, lags, surrogate_sequence, autocorr_bin_width, autocorr_max_lag);
				if (err_ac != SSM_SUCCESS) {
					std::cerr << "ERROR (in spiSeMe_surrogate_sa): error in autocorrelation evaluation.\n";
					return SSM_RUNTIME_ERR;
				}
				permutation_cost = compute_cost_function(autocorrelation_surrogate, autocorrelation_original, cost_code);

				// --- Metropolis step
				delta_E = permutation_cost - previous_cost;
				if ((delta_E > 0) && (gsl_rng_uniform(r) > exp(-delta_E/T))) {
					swap_temp = surrogate_sequence[idx_a];
					surrogate_sequence[idx_a] = surrogate_sequence[idx_b];
					surrogate_sequence[idx_b] = swap_temp;
				} else {
					previous_cost = permutation_cost;
					Ns = Ns + 1;
				}
				Nt = Nt + 1;

				if ((previous_cost <= target_cost) || (Nt == n_total) || (Ns == n_successful))
					break;
			}
			n_swaps = n_swaps + Nt;

			// --- Provide some information
			if (verbose)
				std::cerr << "Temperature: " << T << "\tAccepted swaps: " << Ns << ", Total swaps: " << Nt << "\tCost: " << previous_cost << "\n";

			// --- Cool by cooling factor
			T *= cooling_factor;
			if (previous_cost <= target_cost)
				break;
		}

		if (verbose)
			std::cerr << "Process converged after " << n_swaps << " iterations.\n\n";

		// --- Append sequence to the set of those already generated
		iei_surrogates.push_back(surrogate_sequence);

	}

	return SSM_SUCCESS;
}


double compute_cost_function(const std::vector<double> & autocorr_a, const std::vector<double> & autocorr_b, unsigned int cost_code)
{
	unsigned int K = autocorr_a.size();
	double cost = 0.0;
	if (cost_code == 0) {
		for (int i = 0; i < K; i++) {
			if (fabs(autocorr_a[i] - autocorr_b[i]) > cost) {
			 	cost = fabs(autocorr_a[i] - autocorr_b[i]);
			}
		}
	} else if (cost_code == 1) {
		for (int i = 0; i < K; i++) {
			cost += fabs(autocorr_a[i] - autocorr_b[i]);
		}
	} else if (cost_code == 2) {
		for (int i = 0; i < K; i++) {
			cost += fabs(autocorr_a[i] - autocorr_b[i]) * fabs(autocorr_a[i] - autocorr_b[i]);
		}
	}

	return cost;
}
