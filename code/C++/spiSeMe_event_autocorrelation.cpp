// spiSeMe_event_autocorrelation
//
//	Computes the sample autocorrelation of ieiSequence.
//	The autocorrelation assessment is carried out by binning the lag axis in
//	bins of width bin_width and up to a maximum lag value max_lag.
//	The returned vector "A" contains the counts corresponding to each
//	bin, normalized by N*N*bin_width/T, where N is the number of IEIs
//	and T the sequence duration. The central lag value of each bin is
//	returned in the "lags" vector.
//
//	Both the bin_width (bin width) and max_lag (maximum lag value for
//	which autocorrelation has to be assessed) have to be positive real
//	numbers. In order to provide a meaningful assessment, bin_width must
//	be smaller than max_lag.
//	If max_lag is not an integer multiple of bin_width, the actual
//	maximum lag considered in the assessment of autocorrelation is
//	bin_width*ceil(max_lag/bin_width).
//
//	This function is part of the SpiSeMe package.
//
//
//	REFERENCE:
//	Details on the autocorrelation assessment, including the rationale
//	behind the normalization coefficient, can be found in
//	Chaos 29 (2019), 121102, doi:10.1063/1.5138250
//
//	Citing this source is highly appreciated. Thank you

#ifndef SPISEME_HEADER
	#define SPISEME_HEADER
	#include "spiSeMe.h"
#endif

spiSeMe_return_code spiSeMe_event_autocorrelation(std::vector<double> & A, std::vector<double> & lags, const std::vector<double> & iei_sequence, double bin_width, double max_lag)
{
	// --- Input parsing & validation
	if (bin_width <= 0) {
	 	std::cerr << "ERROR (in spiSeMe_event_autocorrelation): argument 'bin_width' must be positive.\n";
		return SSM_BAD_ARGUMENT;
	}
	if (max_lag <= 0) {
	 	std::cerr << "ERROR (in spiSeMe_event_autocorrelation): argument 'max_lag' must be positive.\n";
		return SSM_BAD_ARGUMENT;
	}
	if (max_lag <= bin_width) {
	 	std::cerr << "ERROR (in spiSeMe_event_autocorrelation): invalid parameters. Bin width must be smaller than maximum lag.\n";
		return SSM_BAD_ARGUMENT;
	}
	unsigned int N = iei_sequence.size();
	if (N < 2) {
		std::cerr << "ERROR (in spiSeMe_event_autocorrelation): input sequences is not an actual vector (length < 2).\n";
		return SSM_BAD_ARGUMENT;
	}
	for (int i = 0; i < N; i++) {
		if (iei_sequence[i] < 0) {
			std::cerr << "ERROR (in spiSeMe_event_autocorrelation): invalid input sequence. One or more IEIs are negative.\n";
			return SSM_BAD_ARGUMENT;
		}
	}

	// --- Compute arrival times out of IEIs
	std::vector<double>	arrival_times(N+1, 0.0);
	arrival_times[0] = 0.0;
	for (int i = 1; i < N + 1; i++) {
		arrival_times[i] = arrival_times[i-1] + iei_sequence[i-1];
	}
	double sequence_duration = arrival_times.back();

	// --- Initialize output arrays (any existing content in A and lags is destroyed)
	unsigned int n_bins = ceil(max_lag / bin_width);
	A.clear();
	lags.clear();
	A.resize(n_bins, 0.0);
	lags.resize(n_bins, 0.0);
	for (int i = 0; i < n_bins; i++) {
		lags[i] = bin_width * (i + 0.5);
	}

	// --- Compute autocorrelation
	double delta_t;
	unsigned int idx;
	for (int i = 1; i < N + 1; i++) {
		for (int j = 0; j < i; j++) {
			delta_t = arrival_times[i] - arrival_times[j];
			// The maximum lag considered is n_bins*bin_width, which is >= max_lag (see initializations)
			if (delta_t < n_bins * bin_width) {
				idx = floor(delta_t / bin_width);
				A[idx]++;
			}
		}
	}

	// --- Normalize autocorrelation [see Chaos 29, 121102 (2019)]
	for (int i = 0; i < n_bins; i++) {
		A[i] /= (N * N * bin_width / sequence_duration);
	}

	return SSM_SUCCESS;
}
