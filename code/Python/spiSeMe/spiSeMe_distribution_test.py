import numpy

def spiSeMe_distribution_test(iei_sequence_A, iei_sequence_B):
	"""
	p, d = spiSeMe_distribution_test(iei_sequence_A, iei_sequence_B)

		Carries out a Kolmogorov-Smirnov test to assess the compatibility between
		the IEI distribution associated to iei_sequence_A and iei_sequence_B.
		The two sequences do not need to have the same number of elements.
		The returned p-value 'p' corresponds to the null-hypothesis of the two
		sample distributions of IEI having the same parent distribution.
		The K-S statistic 'd' is also returned.

	This function is part of the SpiSeMe package.

	REFERENCE:
	The Kolmogorov-Smirnov test and its present implementation are thoroughly
	described in W.H. Press, S.A. Teukolsky W.T. Vetterling and B.P. Flannery,
	Numerical Recipes. The Art of Scientific Computing, 3rd Edition, 2007,
	ISBN 0-521-88068-8.
	"""


	# --- Input validation
	if (not (isinstance(iei_sequence_A,numpy.ndarray)) or not (isinstance(iei_sequence_B,numpy.ndarray))):
			print('ERROR (in spiSeMe_distribution_test): function arguments "iei_sequence_A" and "iei_sequence_B" must be one-dimensional numpy arrays.')
			return False, False
	if ((iei_sequence_A.ndim != 1) or (iei_sequence_B.ndim != 1)):
			print('ERROR (in spiSeMe_distribution_test): function arguments "iei_sequence_A" and "iei_sequence_B" must be one-dimensional numpy arrays.')
			return False, False
	if (numpy.any(iei_sequence_A < 0) or numpy.any(iei_sequence_B < 0)):
		print('ERROR (in spiSeMe_distribution_test): invalid input sequences. One or more IEIs are negative.')
		return False, False

	# --- Key parameters
	n_A = iei_sequence_A.size
	n_B = iei_sequence_B.size
	n_e = numpy.sqrt(n_A * n_B / (n_A + n_B));

	# --- K-S statistic: implementation following Numerical Recipes
	data_A = numpy.sort(iei_sequence_A)
	data_B = numpy.sort(iei_sequence_B)
	cdf_A = 0.0;
	cdf_B = 0.0;
	j1 = 0;
	j2 = 0;
	d = 0;
	while (j1 < n_A and j2 < n_B):
		d1 = data_A[j1]
		d2 = data_B[j2]
		if (d1 <= d2):
			while (j1 < n_A and d1 == data_A[j1]):
				j1 = j1 + 1
				cdf_A = float(j1) / float(n_A);
		if (d2 <= d1):
			while (j2 < n_B and d2 == data_B[j2]):
				j2 = j2 + 1
				cdf_B = float(j2) / float(n_B)
		if (numpy.fabs(cdf_B - cdf_A) > d):
			d = numpy.fabs(cdf_B - cdf_A)

	KS_statistic = (n_e + 0.12 + 0.11/n_e) * d

	def qks(z):
		if (z == 0.):
			q = 1.0
		elif (z < 1.18):
			y = numpy.exp(-1.23370055013616983 / z**2)
			q = 1.0 - 2.25675833419102515 * numpy.sqrt(-numpy.log(y)) * (y + y**9 + y**25 + y**49)
		else:
			y = numpy.exp(-2.0 * z**2)
			q = 2.0 * (y - y**4 + y**9)
		return q

	p = qks(KS_statistic)

	return p, d
