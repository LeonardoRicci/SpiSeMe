// SPISEME_SURROGATE_JODI
//
//	Generates, by means of the JOint DIstribution (JODI) algorithm,
//	surrogate IEI sequences corresponding to the original
//	sequence iei_sequence.
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
//	The JOint DIstribution method (JODI) to generate surrogate event sequences
//	is described in L. Ricci, M. Castelluzzo, L. Minati, and A. Perinelli,
//	Generation of surrogate event sequences via joint distribution of successive
//	inter-event intervals, Chaos 29 (2019), 121102, doi:10.1063/1.5138250
//
//	Citing this source is highly appreciated. Thank you.

#ifndef SPISEME_HEADER
	#define SPISEME_HEADER
	#include "spiSeMe.h"
#endif
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_rng.h>

struct triplet_index_value_cdf {
	unsigned int	index;
	double		value;
	int		rank;
};

bool sort_triplets_by_value(triplet_index_value_cdf a, triplet_index_value_cdf b) { return (a.value < b.value); }
bool sort_triplets_by_index(triplet_index_value_cdf a, triplet_index_value_cdf b) { return (a.index < b.index); }
void build_histogram(gsl_histogram2d*, int, const std::vector<double> &, double, double);
void retrieve_random_bins(int &, int &, gsl_histogram2d*, gsl_rng*, unsigned int, unsigned int);
void retrieve_conditioned_random_bin(int &, gsl_histogram2d*, int, gsl_rng*, unsigned int);

spiSeMe_return_code spiSeMe_surrogate_jodi(std::vector < std::vector<double> > & iei_surrogates, const std::vector<double> & iei_sequence, unsigned int M, bool verbose)
{
	// --- Input parsing & validation
	unsigned int N = iei_sequence.size();
	if (N < 2) {
		std::cerr << "ERROR (in spiSeMe_surrogate_jodi): input iei_sequence is not an actual vector (size < 2).\n";
		return SSM_BAD_ARGUMENT;
	}
	if (M < 1) {
		std::cerr << "ERROR (in spiSeMe_surrogate_jodi): function argument M must be a positive integer.\n";
		return SSM_BAD_ARGUMENT;
	}

	// --- Provide some information
	if (verbose)
		std::cerr << "\n### Starting JODI routine ###\n\n";

	// --- Extract sequence of ranks out of IEI sequence
	std::vector <triplet_index_value_cdf>	triplet_sequence;
	triplet_index_value_cdf			temp_triplet;
	for (int i = 0; i < N; i++) {
		temp_triplet.value = iei_sequence[i];
		temp_triplet.index = i;
		temp_triplet.rank = 0.0;
		triplet_sequence.push_back(temp_triplet);
	}

	std::vector<double>		ranks_sequence;
	sort(triplet_sequence.begin(), triplet_sequence.end(), sort_triplets_by_value);
	for (int i = 0; i < N; i++)
		triplet_sequence[i].rank = i;

	sort(triplet_sequence.begin(), triplet_sequence.end(), sort_triplets_by_index);
	for (int i = 0; i < N; i++)
		ranks_sequence.push_back(triplet_sequence[i].rank);

	// --- Build histograms according to the Freedman-Diaconis rule
	double	isi_min = 0, isi_max = N-1;
	double	bin_size = 0.0 - N/4;
	bin_size += (3*N)/4;
	bin_size *= 2;
	bin_size /= exp(log(N)/3.0);
	unsigned int n_bins = ceil((isi_max - isi_min)/bin_size);

	gsl_histogram2d * histogram_2d = gsl_histogram2d_alloc(n_bins, n_bins);
	gsl_histogram2d_set_ranges_uniform(histogram_2d, isi_min, isi_max, isi_min, isi_max);
	build_histogram(histogram_2d, n_bins, ranks_sequence, isi_min, isi_max);

 	//--- Initialize output array (each of the M generated surrogates is a vector)
	iei_surrogates.clear();
	std::vector<double> temp_iei_surrogate;
	// --- Initialize random engine
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, time(NULL));

 	//--- Generate surrogates: iterate M times
	for (int iter = 0; iter < M; iter++) {
		if (verbose)
			std::cerr << "# Surrogate number " << iter+1 << " out of " << M << ".\n";

		// --- Initialize surrogate ranks generation
		std::vector<double> ranks_surrogate(N, 0.0);
		int	b_xn, b_xnn;
		double	xn_lower, xn_upper, xnn_lower, xnn_upper;

 		//--- Generate the first pair r_1, r_2 according to the sample joint distribution
		retrieve_random_bins(b_xn, b_xnn, histogram_2d, r, n_bins, N);
		gsl_histogram2d_get_xrange(histogram_2d, b_xn, &xn_lower, &xn_upper);
		gsl_histogram2d_get_yrange(histogram_2d, b_xnn, &xnn_lower, &xnn_upper);
		ranks_surrogate[0] = xn_lower + (xn_upper - xn_lower) * gsl_rng_uniform(r);
		ranks_surrogate[1] = xnn_lower + (xnn_upper - xnn_lower) * gsl_rng_uniform(r);

		// --- Iterate the Markov chain to generate successive intervals
		size_t	idx_xn, idx_xnn;
		for (int i = 2; i < N; i++) {
			// --- Compute the conditional distribution P(r_{n+1} | r_n) corresponding to r_n being the last generated element
			gsl_histogram2d_find(histogram_2d, ranks_surrogate[i-1], ranks_surrogate[i-1], &idx_xn, &idx_xnn);
			retrieve_conditioned_random_bin(b_xnn, histogram_2d, idx_xn, r, n_bins);

			// --- Generate the new element according to the conditional distribution P(r_{n+1} | r_n)
			gsl_histogram2d_get_yrange(histogram_2d, b_xnn, &xnn_lower, &xnn_upper);
			ranks_surrogate[i] = xnn_lower + (xnn_upper - xnn_lower) * gsl_rng_uniform(r);
		}

		// --- Transform ranks back into IEIs
		std::vector <triplet_index_value_cdf> triplet_surrogate;
		for (int i = 0; i < N; i++) {
			temp_triplet.index = i;
			temp_triplet.value = ranks_surrogate[i];
			triplet_surrogate.push_back(temp_triplet);
		}
		sort(triplet_sequence.begin(), triplet_sequence.end(), sort_triplets_by_value);
		sort(triplet_surrogate.begin(), triplet_surrogate.end(), sort_triplets_by_value);
		for (int i = 0; i < N; i++)
			triplet_surrogate[i].value = triplet_sequence[i].value;
		sort(triplet_surrogate.begin(), triplet_surrogate.end(), sort_triplets_by_index);

		if (verbose)
			std::cerr << "Process ended; " << N << " elements generated.\n\n";

		// --- Append sequence to the set of those already generated
		temp_iei_surrogate.clear();
		for (int i = 0; i < N; i++)
			temp_iei_surrogate.push_back(triplet_surrogate[i].value);

		iei_surrogates.push_back(temp_iei_surrogate);
	}

	gsl_histogram2d_free(histogram_2d);
	gsl_rng_free(r);

	return SSM_SUCCESS;
}

void build_histogram(gsl_histogram2d * histogram_2d, int n_bins, const std::vector<double> & sequence, double isi_min, double isi_max)
{
	for (int i = 0; i < sequence.size() - 1; i++) {
		if (sequence[i+1] == isi_max)
			gsl_histogram2d_increment(histogram_2d, sequence[i], isi_max - (isi_max - isi_min)/(100.0*n_bins));
		else if (sequence[i] == isi_max)
			gsl_histogram2d_increment(histogram_2d, isi_max - (isi_max - isi_min)/(100.0*n_bins), sequence[i+1]);
		else
			gsl_histogram2d_increment(histogram_2d, sequence[i], sequence[i+1]);
	}

	return;
}

void retrieve_random_bins(int &b_xn, int &b_xnn, gsl_histogram2d * histogram_2d, gsl_rng* r, unsigned int n_bins, unsigned int N)
{
	bool	stop_counting = false;
	int	count_up_to = gsl_rng_uniform_int(r, N); // 0 to N-1
	int	n = 0;

	for (int i = 0; i < n_bins; i++) {
		for (int j = 0; j < n_bins; j++) {
			n += gsl_histogram2d_get(histogram_2d, i, j);
			if (n > count_up_to) {
				b_xn = i;
				b_xnn = j;
				stop_counting = true;
				break;
			}
		}
		if (stop_counting)
			break;
	}

	return;
}

void retrieve_conditioned_random_bin(int &b_xnn, gsl_histogram2d * histogram_2d, int b_xn, gsl_rng* r, unsigned int n_bins)
{
	int	n = 0;

	for (int j = 0; j < n_bins; j++)
		n += gsl_histogram2d_get(histogram_2d, b_xn, j);

	int count_up_to = gsl_rng_uniform_int(r, n); // 0 to n-1
	n = 0;
	for (int j = 0; j < n_bins; j++) {
		n += gsl_histogram2d_get(histogram_2d, b_xn, j);
		if (n > count_up_to) {
			b_xnn = j;
			break;
		}
	}

	return;
}
