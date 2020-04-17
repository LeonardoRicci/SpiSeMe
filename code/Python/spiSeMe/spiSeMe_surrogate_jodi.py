import numpy
import numpy.random

def spiSeMe_surrogate_jodi(iei_sequence, M = 1, verbose = True):
	"""
	iei_surrogates = spiSeMe_surrogate_jodi(iei_sequence, M = 1, verbose = True)

		Generates, by means of the JOint DIstribution (JODI) algorithm,
		surrogate IEI sequences corresponding to the original
		sequence ieiSequence.

	Options:

	'M'		Number of surrogate sequences to generate.
			By default, M = 1. Each of the M surrogate sequences
			is a column in the returned array, which therefore
			has a size L * M, where L is the original sequence
			length.

	'verbose'	Sets the verbosity of the function. If true (default),
			all messages are printed on the command line. If
			false, only critical errors are displayed.

	This function is part of the SpiSeMe package.


	REFERENCE:
	The JOint DIstribution method (JODI) to generate surrogate event sequences
	is described in L. Ricci, M. Castelluzzo, L. Minati, and A. Perinelli,
	Generation of surrogate event sequences via joint distribution of successive
	inter-event intervals, Chaos 29 (2019), 121102, doi:10.1063/1.5138250

	Citing this source is highly appreciated. Thank you.
	"""

	# --- Input validation
	if (not (isinstance(iei_sequence,numpy.ndarray))):
		print('ERROR (in spiSeMe_surrogate_jodi): function argument "iei_sequence" must be a one-dimensional numpy array.')
		return False
	if (iei_sequence.ndim != 1):
		print('ERROR (in spiSeMe_surrogate_jodi): function argument "iei_sequence" must be a one-dimensional numpy array.')
		return False
	L = iei_sequence.size
	try:
		M = int(M)
	except ValueError:
		print('ERROR (in spiSeMe_surrogate_jodi): function argument "M" must be a positive integer or convertible to a positive integer.')
		return False
	if M < 1:
		print('ERROR (in spiSeMe_surrogate_jodi): function argument "M" must be a positive integer or convertible to a positive integer.')
		return False
	if (not (isinstance(verbose,(bool)))):
		print('ERROR (in spiSeMe_surrogate_jodi): function option "verbose" must be Boolean (True / False).')
		return False

	# --- Provide some information
	if verbose is True:
		print('\n### Starting JODI routine ###\n')

	# --- Extract sequence of ranks out of IEI sequence
	order_iei = numpy.argsort(iei_sequence)
	sorted_iei = iei_sequence[order_iei]
	ranks = numpy.linspace(1, L, L)
	ranks_original = numpy.zeros_like(ranks)
	ranks_original[order_iei] = ranks

	# --- Build histograms according to the Freedman-Diaconis rule
	h1counts, h1edges = numpy.histogram(ranks_original, bins='fd')
	nbins = h1counts.size
	h2counts, h2Xedges, h2Yedges = numpy.histogram2d(ranks_original[0:-1], ranks_original[1:], nbins)

	# --- Initialize output array (each of the M generated surrogates is a row).
	iei_surrogates = numpy.empty_like(iei_sequence)

	# --- Generate surrogates: iterate M times
	for iteration_number in range(0, M):
		if verbose is True:
			print('# Surrogate number ', iteration_number + 1, ' out of ', M, '.')

		# --- Initialize surrogate ranks generation
		ranks_surrogate = numpy.zeros_like(ranks_original)
		numpy.random.seed()

		# --- Generate the first pair r_1, r_2 according to the sample joint distribution
		stop_count = numpy.random.randint(1, L+1)	# +1 because upper limit of randint is excluded
		temp_count = 0
		stop_search = False
		for i in range(0, nbins):
			for j in range(0, nbins):
				temp_count = temp_count + h2counts[i, j]
				if (temp_count >= stop_count):
					stop_search = True
					break
			if (stop_search):
				break
		ranks_surrogate[0] = h2Xedges[i] + numpy.random.rand()*(h2Xedges[i+1] - h2Xedges[i])
		ranks_surrogate[1] = h2Yedges[j] + numpy.random.rand()*(h2Yedges[j+1] - h2Yedges[j])

		# --- Iterate the Markov chain to generate successive intervals
		n = 2
		while n < L:
			# --- Compute the conditional distribution P(r_n | r_{n-1}) corresponding to r_{n-1} being the last generated element
			last_extracted_bin = numpy.digitize(ranks_surrogate[n-1], h2Xedges) - 1	# -1 because numpy.digitize() starts counting bins from 1
			conditional_distribution = h2counts[last_extracted_bin]
			conditional_distribution_normalization = numpy.sum(conditional_distribution)

			# --- Generate the new element according to the conditional distribution P(r_n | r_{n-1})
			stop_count = numpy.random.randint(1, conditional_distribution_normalization+1) # +1 because upper limit of randint is excluded
			temp_count = 0
			for j in range(0, nbins):
				temp_count = temp_count + conditional_distribution[j]
				if (temp_count >= stop_count):
					break
			ranks_surrogate[n] = h2Yedges[j] + numpy.random.rand()*(h2Yedges[j+1] - h2Yedges[j])
			n = n + 1

		# --- Transform ranks back into IEIs
		order_ranks = numpy.argsort(ranks_surrogate)
		iei_surrogate_seq = numpy.zeros_like(ranks_surrogate)
		iei_surrogate_seq[order_ranks] = sorted_iei

		if verbose is True:
			print('Process ended; ', n, ' elements generated.\n')

		# --- Append sequence to the set of those already generated
		if (iteration_number == 0):
			iei_surrogates = iei_surrogate_seq
		else:
			iei_surrogates = numpy.vstack((iei_surrogates, iei_surrogate_seq))

	return iei_surrogates
