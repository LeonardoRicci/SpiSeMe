import numpy
import numpy.random
import spiSeMe.spiSeMe_event_autocorrelation as spiseme_ac

def spiSeMe_surrogate_sa(iei_sequence, autocorr_bin_width, autocorr_max_lag, a = 0.9, T = 0.1, C = -1, cost_function = 'max', n_total = -1, n_successful = -1, M = 1, verbose = True):
	"""
	iei_surrogates = spiSeMe_surrogate_sa(iei_sequence, autocorr_bin_width, autocorr_max_lag, a = 0.9, T = 0.1, C = -1, cost_function = 'max', n_total = -1, n_successful = -1, M = 1, verbose = True)

		Generates, by means of the Simulated Annealing (SA) algorithm,
		surrogate IEI sequences corresponding to the original
		sequence ieiSequence.

		The evaluation of the cost function relies on the assessment of
		event autocorrelation. This assessment is implemented by the
		function spiSeMe_event_autocorrelation within this package. The mandatory
		parameters 'autocorr_bin_width', 'autocorr_max_lag' correspond to
		the 'bin_width', 'max_lag' parameters of the spiSeMe_event_autocorrelation
		function: please refer to its documentation for further details.

	Options:

	'a'		Cooling factor (in references, alpha). Default 0.9.

	'T'		Starting temperature. Default 0.1.

	'C'		Target cost below which the routine is stopped. If not
			specified, the function attempts to estimate a target cost as
			the standard deviation of the original sequence's
			autocorrelation. This estimate is not necessarily a good
			choice for every sequence. Similarly, different metrics (see
			'costFunction') will lead to a different accuracy under the
			default target cost.

	'cost_function'	Sets the metric used to evaluate the cost function:
			'max': Cost equals the maximim absolute difference between
			original and surrogate autocorrelations (default setting).
			'L1': Cost equals the sum over all bins of the absolute
			differences between original and surrogate autocorrelations.
			'L2': Cost equals the sum over all bins of the squared
			differences between original and surrogate autocorrelations.

	'n_total' Number of swaps to be carried out before cooling.
			Default is L, where L is the sequence length.

	'n_successful' Number of accepted swaps to be carried out before cooling.
			Default is L/2, where L is the sequence length.

	'M'		Number of surrogate sequences to be generated. By
			default, M = 1. Each of the M surrogate sequences
			is a column in the returned array, which therefore
			has size MxL, where L is the original sequence
			length.

	'verbose' Sets the verbosity of the function. If true (default),
			all messages are displayed. If false, only critical errors
			are displayed.

	This function is part of the SpiSeMe package.


	REFERENCE:
	The simulated annealing algorithm for surrogate generation was originally
	proposed by T. Schreiber in Phys. Rev. Lett. 80 (1998), 2105,
	doi:10.1103/PhysRevLett.80.2105

	Details on the autocorrelation assessment, as well as a comparative study
	of surrogate generation algorithms, can be found in
	Chaos 29 (2019), 121102, doi:10.1063/1.5138250

	Citing this source is highly appreciated. Thank you.

	See also spiSeMe_event_autocorrelation
	"""

	# --- Input validation
	if (not (isinstance(iei_sequence,numpy.ndarray))):
			print('ERROR (in spiSeMe_surrogate_sa): function argument "iei_sequence" must be a one-dimensional numpy array.')
			return False
	if (iei_sequence.ndim != 1):
			print('ERROR (in spiSeMe_surrogate_sa): function argument "iei_sequence" must be a one-dimensional numpy array.')
			return False
	L = iei_sequence.size
	try:
		M = int(M)
	except ValueError:
		print('ERROR (in spiSeMe_surrogate_sa): function argument "M" must be a positive integer or convertible to a positive integer.')
		return False
	if M < 1:
		print('ERROR (in spiSeMe_surrogate_sa): function argument "M" must be a positive integer or convertible to a positive integer.')
		return False
	if (not (isinstance(verbose,(bool)))):
		print('ERROR (in spiSeMe_surrogate_sa): function option "verbose" must be a boolean.')
		return False
	try:
		autocorr_bin_width = float(autocorr_bin_width)
	except:
		print('ERROR (in spiSeMe_surrogate_sa): function argument "autocorr_bin_width" must be a float or convertible to a float.')
		return False
	try:
		autocorr_max_lag = float(autocorr_max_lag)
	except:
		print('ERROR (in spiSeMe_surrogate_sa): function argument "autocorr_max_lag" must be a float or convertible to a float.')
		return False
	if (autocorr_max_lag <= 0) or (autocorr_bin_width <= 0):
		print('ERROR (in spiSeMe_surrogate_sa): function arguments "autocorr_bin_width" and "max_lag" must be positive.')
		return False
	if (autocorr_max_lag <= autocorr_bin_width):
		print('ERROR (in spiSeMe_surrogate_sa): invalid parameter. "autocorr_bin_width" must be smaller than "autocorr_max_lag".')
		return False
	try:
		starting_temperature = float(T)
	except ValueError:
		print('ERROR (in spiSeMe_surrogate_sa): function argument "starting_temperature" must be a float or convertible to a float.')
		return False
	try:
		cooling_factor = float(a)
	except ValueError:
		print('ERROR (in spiSeMe_surrogate_sa): function argument "cooling_factor" must be a float or convertible to a float.')
		return False
	if (starting_temperature < 0):
		print('ERROR (in spiSeMe_surrogate_sa): invalid parameter. "starting_temperature" must be positive.')
		return False
	if ((cooling_factor >= 1.0) or (cooling_factor <= 0)):
		print('ERROR (in spiSeMe_surrogate_sa): invalid parameter. "cooling_factor" must be positive and less than unity.')
		return False
	try:
		n_total = int(n_total)
	except ValueError:
		print('ERROR (in spiSeMe_surrogate_sa): function argument "n_total" must be an integer.')
		return False
	if (n_total <= 0):
		n_total = int(L)
	try:
		n_successful = int(n_successful)
	except ValueError:
		print('ERROR (in spiSeMe_surrogate_sa): function argument "n_successful" must be an integer.')
		return False
	if (n_successful <= 0):
		n_successful = int(L / 2)

	if ((cost_function != 'max') and (cost_function != 'L1') and (cost_function != 'L2')):
		print('ERROR (in spiSeMe_surrogate_sa): function argument "cost_function" is not among the possible ones ("max", "L1" or "L2").')
		return False

	# --- Assess autocorrelation of original sequence
	original_autocorr, lags = spiseme_ac.spiSeMe_event_autocorrelation(iei_sequence, autocorr_bin_width, autocorr_max_lag)

	# --- If not assigned, estimate the target cost as the standard deviation of the Autocorrelation divided by the square root of its length
	if (C <= 0):
		target_cost = numpy.std(original_autocorr)
		print('WARNING: Automatic setting of target cost might be too small for the routine to converge in a reasonable time.')
	else:
		target_cost = C

	# --- Provide some information
	if verbose is True:
		print('\n### Starting simulated annealing routine ###')
		print('###\tTarget cost: ', target_cost)
		print('###\tCooling factor: ', cooling_factor)
		print('###\tStarting T: ', starting_temperature)
		print('###\tSuccessful swaps before cooling: ', n_successful)
		print('###\tTotal swaps before cooling: ', n_total, '\n')

	# --- Initialize output array (each of the M generated surrogates is a row)
	iei_surrogates = numpy.empty_like(iei_sequence)

	# --- Cost function definition
	def compute_cost_function_max(autocorr_surr, autocorr_original):
		distance = numpy.amax(numpy.fabs(autocorr_surr - autocorr_original))
		return distance
	def compute_cost_function_L1(autocorr_surr, autocorr_original):
		distance = numpy.sum(numpy.fabs(autocorr_surr - autocorr_original))
		return distance
	def compute_cost_function_L2(autocorr_surr, autocorr_original):
		distance = numpy.sum(numpy.square(autocorr_surr - autocorr_original))
		return distance

	# --- Generate surrogates: iterate M times
	for iter in range(0, M):

		if verbose is True:
			print('# Surrogate number ', iter + 1, ' out of ', M, '.')

		# --- Starting conditions
		numpy.random.seed()
		T = starting_temperature
		run_cooling = True
		run_randomization = True
		n_swaps = 0
		sequence_surrogate_IEI = numpy.copy(iei_sequence)
		numpy.random.shuffle(sequence_surrogate_IEI)
		surrogate_autocorr, lags = spiseme_ac.spiSeMe_event_autocorrelation(sequence_surrogate_IEI, autocorr_bin_width, autocorr_max_lag);

		if (cost_function == 'max'):
			previous_cost = compute_cost_function_max(surrogate_autocorr, original_autocorr)
		elif (cost_function == 'L1'):
			previous_cost = compute_cost_function_L1(surrogate_autocorr, original_autocorr)
		elif (cost_function == 'L2'):
			previous_cost = compute_cost_function_L2(surrogate_autocorr, original_autocorr)

		# --- Outermost iteration: cooling stages
		while (True):

			Ns = 0;
			Nt = 0;

			# --- Innermost iteration: swaps & Metropolis step
			while (True):

				idxA = numpy.random.randint(0, high=L)
				idxB = numpy.random.randint(0, high=L)

				sequence_surrogate_IEI[[idxA, idxB]] = sequence_surrogate_IEI[[idxB, idxA]]

				surrogate_autocorr, lags = spiseme_ac.spiSeMe_event_autocorrelation(sequence_surrogate_IEI, autocorr_bin_width, autocorr_max_lag)
				if (cost_function == 'max'):
					permutation_cost = compute_cost_function_max(surrogate_autocorr, original_autocorr)
				elif (cost_function == 'L1'):
					permutation_cost = compute_cost_function_L1(surrogate_autocorr, original_autocorr)
				elif (cost_function == 'L2'):
					permutation_cost = compute_cost_function_L2(surrogate_autocorr, original_autocorr)

				# --- Metropolis step
				deltaE = permutation_cost - previous_cost
				if ((deltaE > 0) and (numpy.random.rand(1) > numpy.exp(- deltaE / T))):
					sequence_surrogate_IEI[[idxA, idxB]] = sequence_surrogate_IEI[[idxB, idxA]]
				else:
					previous_cost = permutation_cost
					Ns += 1
				Nt += 1

				if ((previous_cost <= target_cost) or (Nt == n_total) or (Ns == n_successful)):
					break

			n_swaps = n_swaps + Nt;

			# --- Provide some information
			if verbose is True:
				print('Temperature: ', T, '\tAccepted swaps: ', Ns, ', Total swaps: ', Nt, '\tCost: ', previous_cost)

			# --- Cool by cooling factor
			T = cooling_factor * T;
			if (previous_cost <= target_cost):
				break

		if verbose is True:
			print('Target cost reached after ', n_swaps, ' total swaps.\n')

		# --- Append sequence to the set of those already generated
		if (iter == 0):
			iei_surrogates = sequence_surrogate_IEI
		else:
			iei_surrogates = numpy.vstack((iei_surrogates, sequence_surrogate_IEI))

	return iei_surrogates
