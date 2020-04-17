import numpy

def spiSeMe_event_cross_correlation(iei_sequence_A, iei_sequence_B, bin_width, max_lag):
	"""
	C, lags = spiSeMe_event_cross_correlation(iei_sequence_A, iei_sequence_B, bin_width, max_lag)

		Computes the sample cross-correlation between iei_sequence_A and iei_sequence_B.

		The cross correlation assessment is carried out by binning the lag axis in
		bins of width bin_width for lag values from -max_lag to max_lag.
		The returned vector "C" contains the counts corresponding to each
		bin. The central lag value of each bin is returned in the "lags" vector.

		Both the bin_width (bin width) and max_lag (maximum lag value for
		which autocorrelation has to be assessed) have to be positive real
		numbers. In order to provide a meaningful assessment, bin_width must
		be smaller than max_lag.
		If max_lag is not an half-integer multiple of bin_width, the actual
		maximum lag considered in the assessment of autocorrelation is
		bin_width*(ceil(max_lag/bin_width) + 0.5).

	This function is part of the SpiSeMe package.

	"""


	# --- Input validation
	if ((not (isinstance(iei_sequence_A,numpy.ndarray))) or (not (isinstance(iei_sequence_B,numpy.ndarray)))):
			print('ERROR (in spiSeMe_event_cross_correlation): function arguments "iei_sequence_A" and "iei_sequence_B" must be one-dimensional numpy arrays.')
			return False, False
	if ((iei_sequence_A.ndim != 1) or (iei_sequence_B.ndim != 1)):
			print('ERROR (in spiSeMe_event_cross_correlation): function arguments "iei_sequence_A" and "iei_sequence_B" must be one-dimensional numpy arrays.')
			return False, False
	Na = iei_sequence_A.size
	Nb = iei_sequence_B.size
	try:
		bin_width = float(bin_width)
	except:
		print('ERROR (in spiSeMe_event_cross_correlation): function argument "bin_width" must be a float or convertible to a float.')
		return False, False
	try:
		max_lag = float(max_lag)
	except:
		print('ERROR (in spiSeMe_event_cross_correlation): function argument "max_lag" must be a float or convertible to a float.')
		return False, False
	if (max_lag <= 0) or (bin_width <= 0):
		print('ERROR (in spiSeMe_event_cross_correlation): function arguments "bin_width" and "max_lag" must be positive.')
		return False, False
	if (max_lag <= bin_width):
		print('ERROR (in spiSeMe_event_cross_correlation): invalid parameters. "bin_width" must be smaller than "max_lag".')
		return False, False
	if (numpy.any(iei_sequence_A < 0) or numpy.any(iei_sequence_B < 0)):
		print('ERROR (in spiSeMe_event_cross_correlation): invalid input sequences. One or more IEIs are negative.')
		return False, False

	# --- Compute arrival times out of IEIs
	arrival_times_a = numpy.cumsum(iei_sequence_A)
	arrival_times_a = numpy.insert(arrival_times_a, 0, 0.0)
	arrival_times_b = numpy.cumsum(iei_sequence_B)
	arrival_times_b = numpy.insert(arrival_times_b, 0, 0.0)

	# --- Initialize output arrays
	n_bins = int(2 * numpy.ceil(max_lag / bin_width) + 1)
	C = numpy.zeros(n_bins)
	lags = bin_width * numpy.arange(start = -int(numpy.floor(n_bins/2)), stop = int(numpy.floor(n_bins/2)) + 1)
	actual_max_lag = lags[-1] + 0.5*bin_width;

	# --- Compute autocorrelation
	for i in range(0, Na + 1):
	 	at_temp = arrival_times_a[i] - arrival_times_b
	 	idx = numpy.floor(at_temp[numpy.fabs(at_temp) < actual_max_lag] / bin_width + 0.5) + numpy.floor(n_bins/2)
	 	for j in range(0, idx.size):
	 		C[int(idx[j])] += 1

	return C, lags
