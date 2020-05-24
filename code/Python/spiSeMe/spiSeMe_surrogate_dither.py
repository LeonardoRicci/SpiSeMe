import numpy
import numpy.random

def spiSeMe_surrogate_dither(iei_sequence, dither_distribution = 'uniform', distribution_parameter = -1, M = 1, verbose = True):
	"""
	Generates surrogates of Inter-Event-Intervals (IEI) sequences through dithering.

	iei_surrogates = spiSeMe_surrogate_dither(iei_sequence, dither_distribution = 'uniform', distribution_parameter = -1, M = 1, verbose = True)

		Generates, by means of dithering, surrogate IEI sequences corresponding
		to the original sequence ieiSequence.
		Dithering consists in adding random time intervals to the time
		coordinate of each event in the sequence. Random intervals follow a
		zero-mean distribution whose shape and width can be set by the
		user.

	Options:

	'dither_distribution' Specifies from which distribution the
			dithering r.v. are drawn. Can either be 'uniform' (default),
			'normal', or 'triangular'.

	'distribution_parameter' In the case of "uniform" or "triangular", sets
			the half-width of the distribution from	which dithering
			elements have to be drawn (see 'dither_distribution').
			In the case of "normal", sets the distribution standard
			deviation. By default, D is estimated as half the minimum
			IEI within the original sequence.

	'M'		Number of surrogate sequences to be generated. By
			default, M = 1. Each of the M surrogate sequences
			is a column in the returned array, which therefore
			has size MxL, where L is the original sequence
			length.

	'verbose' Sets the verbosity of the function. If True (default),
			all messages are displayed. If False, only critical errors
			are displayed.

	This function is part of the SpiSeMe package.


	REFERENCE:
	A comparative study of surrogate generation algorithms can be found in
	Chaos 29 (2019), 121102, doi:10.1063/1.5138250

	Citing this source is highly appreciated. Thank you.
	"""

	# --- Input validation
	if (not (isinstance(iei_sequence,numpy.ndarray))):
			print('ERROR (in spiSeMe_surrogate_dither): function argument "iei_sequence" must be a one-dimensional numpy array.')
			return False
	if (iei_sequence.ndim != 1):
			print('ERROR (in spiSeMe_surrogate_dither): function argument "iei_sequence" must be a one-dimensional numpy array.')
			return False
	L = iei_sequence.size
	try:
		M = int(M)
	except ValueError:
		print('ERROR (in spiSeMe_surrogate_dither): function argument "M" must be a positive integer or convertible to a positive integer.')
		return False
	if M < 1:
		print('ERROR (in spiSeMe_surrogate_dither): function argument "M" must be a positive integer or convertible to a positive integer.')
		return False
	if (not (isinstance(verbose,(bool)))):
		print('ERROR (in spiSeMe_surrogate_dither): function option "verbose" must be a boolean.')
		return False
	try:
		distribution_parameter = float(distribution_parameter)
	except:
		print('ERROR (in spiSeMe_surrogate_dither): function argument "distribution_parameter" must be a float or convertible to a float.')
		return False

	if ((dither_distribution != 'uniform') and (dither_distribution != 'normal') and (dither_distribution != 'triangular')):
		print('ERROR (in spiSeMe_surrogate_dither): function argument "dither_distribution" does not match any of the allowed values ("uniform", "normal", "triangular").')
		return False


	# --- If not assigned, estimate the distribution parameter as half the minimum IEI (in this way, in the case of a uniform distribution, we avoid to change the order of events).
	if (distribution_parameter <= 0):
		distribution_parameter = 0.5 * numpy.min(iei_sequence)

	# --- Provide some information
	if verbose is True:
		print('\n### Starting dithering routine ###')
		print('###\tDither distribution:', dither_distribution)
		if (dither_distribution == 'uniform'):
			print('###\t\tU(-D,D), with D =', distribution_parameter)
		elif (dither_distribution == 'normal'):
			print('###\t\tN(0,s^2), with s =', distribution_parameter)
		elif (dither_distribution == 'triangular'):
			print('###\t\tT(-D,D), with D =', distribution_parameter)

	# --- Build arrival times to be dithered
	original_arrival_times = numpy.concatenate((numpy.array([[0]]), numpy.cumsum(iei_sequence)), axis=None)

	# --- Initialize output array (each of the M generated surrogates is a row)
	iei_surrogates = numpy.empty_like(iei_sequence)

	# --- Generate surrogates: iterate M times
	for iter in range(0, M):

		if verbose is True:
			print('# Surrogate number ', iter + 1, ' out of ', M, '.')

		# --- Initialize surrogate ranks generation
		numpy.random.seed()
		arrival_times = numpy.copy(original_arrival_times)

		# --- Apply dithering to arrival times
		if (dither_distribution == 'uniform'):
			deltas = distribution_parameter * (-1.0 + 2.0 * numpy.random.rand(arrival_times.size))
			arrival_times = arrival_times + deltas;
		elif (dither_distribution == 'normal'):
			deltas = distribution_parameter * numpy.random.normal(0.0, 1.0, arrival_times.size)
			arrival_times = arrival_times + deltas;
		elif (dither_distribution == 'triangular'):
			deltas = distribution_parameter * (numpy.random.rand(arrival_times.size) - numpy.random.rand(arrival_times.size))
			arrival_times = arrival_times + deltas;

		# --- Transform ranks back into IEIs
		this_surrogate = numpy.copy(arrival_times[1:arrival_times.size] - arrival_times[0:arrival_times.size-1])

		if verbose is True:
			print('Process ended.\n')

		# --- Append sequence to the set of those already generated
		if (iter == 0):
			iei_surrogates = numpy.copy(this_surrogate)
		else:
			iei_surrogates = numpy.vstack((iei_surrogates, this_surrogate))

	return iei_surrogates
