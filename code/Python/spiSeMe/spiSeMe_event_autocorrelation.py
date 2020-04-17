import numpy

def spiSeMe_event_autocorrelation(iei_sequence, bin_width, max_lag):
	"""
	A, lags = spiSeMe_event_autocorrelation(iei_sequence, bin_width, max_lag)

		Computes the sample autocorrelation of ieiSequence.
		The autocorrelation assessment is carried out by binning the lag axis in
		bins of width bin_width and up to a maximum lag value max_lag.
		The returned array "A" contains the counts corresponding to each
		bin, normalized by N*N*bin_width/T, where N is the number of IEIs
		and T the sequence duration. The central lag value of each bin is
		returned in the "lags" array.

		Both the bin_width (bin width) and max_lag (maximum lag value for
		which autocorrelation has to be assessed) have to be positive real
		numbers. In order to provide a meaningful assessment, bin_width must
		be smaller than max_lag.
		If max_lag is not an integer multiple of bin_width, the actual
		maximum lag considered in the assessment of autocorrelation is
		bin_width*ceil(max_lag/bin_width).

	This function is part of the SpiSeMe package.


	REFERENCE:
	Details on the autocorrelation assessment, including the rationale
	behind the normalization coefficient, can be found in
	Chaos 29 (2019), 121102, doi:10.1063/1.5138250

	Citing this source is highly appreciated. Thank you.
	"""


	# --- Input validation
	if (not (isinstance(iei_sequence,numpy.ndarray))):
			print('ERROR (in spiSeMe_event_autocorrelation): function argument "iei_sequence" must be a one-dimensional numpy array.')
			return False, False
	if (iei_sequence.ndim != 1):
			print('ERROR (in spiSeMe_event_autocorrelation): function argument "iei_sequence" must be a one-dimensional numpy array.')
			return False, False
	N = iei_sequence.size
	try:
		bin_width = float(bin_width)
	except:
		print('ERROR (in spiSeMe_event_autocorrelation): function argument "bin_width" must be a float or convertible to a float.')
		return False, False
	try:
		max_lag = float(max_lag)
	except:
		print('ERROR (in spiSeMe_event_autocorrelation): function argument "max_lag" must be a float or convertible to a float.')
		return False, False
	if (max_lag <= 0) or (bin_width <= 0):
		print('ERROR (in spiSeMe_event_autocorrelation): function arguments "bin_width" and "max_lag" must be positive.')
		return False, False
	if (max_lag <= bin_width):
		print('ERROR (in spiSeMe_event_autocorrelation): invalid parameters. "bin_width" must be smaller than "max_lag".')
		return False, False
	if (numpy.any(iei_sequence < 0)):
		print('ERROR (in spiSeMe_event_autocorrelation): invalid input sequence. One or more IEIs are negative.')
		return False, False

	# --- Compute arrival times out of IEIs
	arrival_times = numpy.cumsum(iei_sequence)
	arrival_times = numpy.insert(arrival_times, 0, 0.0)
	sequence_duration = arrival_times[-1]

	# --- Initialize output arrays
	A = numpy.zeros(int(numpy.floor(max_lag / bin_width)))
	lags = bin_width * (numpy.arange(start = 0, stop = A.size) + 0.5)

	# --- Compute autocorrelation
	for i in range(1, arrival_times.size):
	 	delta_arrival_times = arrival_times[i] - arrival_times[0:i]
	 	idx = numpy.floor(delta_arrival_times[delta_arrival_times < bin_width*A.size] / bin_width)
	 	# The maximum lag considered is nBins*binWidth, which is >= maxLag (see initializations)
	 	for j in range(0, idx.size):
	 		A[int(idx[j])] += 1

	# --- Normalize autocorrelation [see Chaos 29, 121102 (2019)]
	A = A / (N * N * bin_width / sequence_duration)

	return A, lags
