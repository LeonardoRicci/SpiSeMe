// spiSeMe_event_cross_correlation
//
//	Computes the sample cross correlation between iei_sequence_a and iei_sequence_b.
//	The cross correlation assessment is carried out by binning the lag axis in
//	bins of width bin_width for lag values from -max_lag to max_lag.
//	The returned vector "C" contains the counts corresponding to each
//	bin. The central lag value of each bin is returned in the "lags" vector.
//
//	Both the bin_width (bin width) and max_lag (maximum lag value for
//	which autocorrelation has to be assessed) have to be positive real
//	numbers. In order to provide a meaningful assessment, bin_width must
//	be smaller than max_lag.
//	If max_lag is not an half-integer multiple of bin_width, the actual
//	maximum lag considered in the assessment of autocorrelation is
//	bin_width*(ceil(max_lag/bin_width) + 0.5).
//
//	This function is part of the SpiSeMe package.

#ifndef SPISEME_HEADER
	#define SPISEME_HEADER
	#include "spiSeMe.h"
#endif

spiSeMe_return_code spiSeMe_event_cross_correlation(std::vector<double> & C, std::vector<double> & lags, const std::vector<double> & iei_sequence_A, const std::vector<double> & iei_sequence_B, double bin_width, double max_lag)
{
	// --- Input parsing & validation
	if (bin_width <= 0) {
	 	std::cerr << "ERROR (in spiSeMe_event_cross_correlation): argument 'bin_width' must be positive.\n";
		return SSM_BAD_ARGUMENT;
	}
	if (max_lag <= 0) {
	 	std::cerr << "ERROR (in spiSeMe_event_cross_correlation): argument 'max_lag' must be positive.\n";
		return SSM_BAD_ARGUMENT;
	}
	if (max_lag <= bin_width) {
	 	std::cerr << "ERROR (in spiSeMe_event_cross_correlation): invalid parameters. Bin width must be smaller than maximum lag.\n";
		return SSM_BAD_ARGUMENT;
	}
	unsigned int Na = iei_sequence_A.size();
	unsigned int Nb = iei_sequence_B.size();
	unsigned int N = (Na > Nb)? Na : Nb;
	for (int i = 0; i < Na; i++) {
		if (iei_sequence_A[i] < 0) {
			std::cerr << "ERROR (in spiSeMe_event_cross_correlation): invalid input sequence. One or more IEIs are negative.\n";
			return SSM_BAD_ARGUMENT;
		}
	}
	for (int i = 0; i < Nb; i++) {
		if (iei_sequence_B[i] < 0) {
			std::cerr << "ERROR (in spiSeMe_event_cross_correlation): invalid input sequence. One or more IEIs are negative.\n";
			return SSM_BAD_ARGUMENT;
		}
	}

	// --- Compute arrival times out of IEIs
	std::vector<double>	arrival_times_A(Na+1, 0.0);
	std::vector<double>	arrival_times_B(Nb+1, 0.0);
	arrival_times_A[0] = 0.0;
	arrival_times_B[0] = 0.0;
	for (int i = 1; i < N + 1; i++) {
		if (i < Na + 1)
			arrival_times_A[i] = arrival_times_A[i-1] + iei_sequence_A[i-1];
		if (i < Nb + 1)
			arrival_times_B[i] = arrival_times_B[i-1] + iei_sequence_B[i-1];
	}

	// --- Initialize output arrays (any existing content in C and lags is destroyed)
	int n_bins = 2 * ceil(max_lag / bin_width) + 1;
	int n_bins_half = floor(n_bins / 2);
	C.clear();
	lags.clear();
	for (int i = -n_bins_half; i <= n_bins_half; i++) {
		lags.push_back(bin_width * i);
		C.push_back(0.0);
	}
	double actual_max_lag = lags.back() + 0.5*bin_width;

	// --- Compute cross correlation
	double delta_t;
	int idx;
	for (int i = 0; i < Na + 1; i++) {
		for (int j = 0; j < Nb + 1; j++) {
			delta_t = arrival_times_A[i] - arrival_times_B[j];
			if (fabs(delta_t) < actual_max_lag) {
				idx = floor(delta_t/bin_width + 0.5) + n_bins_half;
				C[idx]++;
			}
		}
	}

	return SSM_SUCCESS;
}
