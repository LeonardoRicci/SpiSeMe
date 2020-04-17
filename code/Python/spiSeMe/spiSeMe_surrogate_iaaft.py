import numpy
import numpy.random
import numpy.fft

def spiSeMe_surrogate_iaaft(iei_sequence, exactly_preserve = 'distribution', M = 1, verbose = True):
	"""
	iei_surrogates = spiSeMe_surrogate_iaaft(iei_sequence, exactly_preserve = 'distribution', M = 1, verbose = True)

		Generates, by means of the Iterative Amplitude Adjusted Fourier Transform
		(IAAFT) algorithm, surrogate IEI sequences corresponding to the original
		sequence ieiSequence.

	Options:

	'exactly_preserve'	Specifies which function has to be preserved
			exactly by the surrogate sequence. Can either be
			'distribution' (default) or 'spectrum'.

	'M'		Number of surrogate sequences to be generated. By
			default, M = 1. Each of the M surrogate sequences
			is a row in the returned array, which therefore
			has size MxL, where L is the original sequence
			length.

	'verbose' Sets the verbosity of the function. If True (default),
			all messages are displayed. If False, only critical errors
			are displayed.

	This function is part of the SpiSeMe package.


	REFERENCE:
	The IAAFT (Iterative Amplitude-Adjusted Fourier Transform) algorithm
	for surrogate generation was originally proposed by T. Schreiber
	and A. Schmitz in Phys. Rev. Lett. 77 (1996), 635,
	doi:10.1103/PhysRevLett.77.635

	A comparative study of surrogate generation algorithms can be found in
	Chaos 29 (2019), 121102, doi:10.1063/1.5138250

	Citing this source is highly appreciated. Thank you.
	"""


	# --- Input validation
	if (not (isinstance(iei_sequence,numpy.ndarray))):
			print('ERROR (in spiSeMe_surrogate_iaaft): function argument "iei_sequence" must be a one-dimensional numpy array.')
			return False
	if (iei_sequence.ndim != 1):
			print('ERROR (in spiSeMe_surrogate_iaaft): function argument "iei_sequence" must be a one-dimensional numpy array.')
			return False
	L = iei_sequence.size
	try:
		M = int(M)
	except ValueError:
		print('ERROR (in spiSeMe_surrogate_iaaft): function argument "M" must be a positive integer or convertible to a positive integer.')
		return False
	if M < 1:
		print('ERROR (in spiSeMe_surrogate_iaaft): function argument "M" must be a positive integer or convertible to a positive integer.')
		return False
	if (not (isinstance(verbose,(bool)))):
		print('ERROR (in spiSeMe_surrogate_iaaft): function option "verbose" must be a boolean.')
		return False
	if ((exactly_preserve != 'distribution') and (exactly_preserve != 'spectrum')):
		print('ERROR (in spiSeMe_surrogate_iaaft): function option "exactly_preserve" can only be "distribution" or "spectrum".')
		return False

	# --- Provide some information
	if verbose is True:
		print('\n### Starting IAAFT routine ###')
		print('###\tExactly preserving: ', exactly_preserve, '\n')

	# --- Assess spectrum & distribution of original sequence
	spectrum_magnitude_original = numpy.absolute(numpy.fft.rfft(iei_sequence))
	distribution_original = numpy.sort(iei_sequence)

 	# --- Initialize output array (each of the M generated surrogates is a row)
	iei_surrogates = numpy.empty_like(iei_sequence);

	# --- Define check of convergence (sequence did not change with respect to last iteration)
	def check_convergence(sequence_this, sequence_prev):
		delta = numpy.sum(numpy.fabs(sequence_prev - sequence_this))
		if (delta == 0):
			return True
		else:
			return False

	# --- Generate surrogates: iterate M times
	for iter in range(0, M):

		if verbose is True:
			print('# Surrogate number ', iter + 1, ' out of ', M, '.')

		# --- Starting conditions
		numpy.random.seed()
		run_iterations = True
		sequence_surrogate_IEI = numpy.copy(iei_sequence)
		numpy.random.shuffle(sequence_surrogate_IEI)
		sequence_surrogate_IEI_previous = numpy.copy(sequence_surrogate_IEI)

		# --- Iterative algorithm (switched depending on target function to be matched exactly.
		n_iter = 0
		while (run_iterations == True):

			n_iter = n_iter + 1;

			spectrum_surrogate = numpy.fft.rfft(sequence_surrogate_IEI);
			phases_surrogate = numpy.angle(spectrum_surrogate);
			spectrum_surrogate = spectrum_magnitude_original * (numpy.cos(phases_surrogate) + 1j*numpy.sin(phases_surrogate));
			sequence_surrogate_IEI = numpy.fft.irfft(spectrum_surrogate, L);

			order_elements = numpy.argsort(sequence_surrogate_IEI)
			sequence_surrogate_IEI[order_elements] = distribution_original

			if (check_convergence(sequence_surrogate_IEI, sequence_surrogate_IEI_previous)):
				if (exactly_preserve == 'spectrum'):
					spectrum_surrogate = numpy.fft.rfft(sequence_surrogate_IEI);
					phases_surrogate = numpy.angle(spectrum_surrogate);
					spectrum_surrogate = spectrum_magnitude_original * (numpy.cos(phases_surrogate) + 1j*numpy.sin(phases_surrogate));
					sequence_surrogate_IEI = numpy.fft.irfft(spectrum_surrogate, L);
				break

			sequence_surrogate_IEI_previous = numpy.copy(sequence_surrogate_IEI)

		if verbose is True:
			print('Process converged after ', n_iter, ' iterations.\n')

		# --- Append sequence to the set of those already generated
		if (iter == 0):
			iei_surrogates = sequence_surrogate_IEI
		else:
			iei_surrogates = numpy.vstack((iei_surrogates, sequence_surrogate_IEI))

	return iei_surrogates
