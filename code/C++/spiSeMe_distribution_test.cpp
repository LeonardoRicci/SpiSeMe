// SPISEME_DISTRIBUTION_TEST
//
//	Carries out a Kolmogorov-Smirnov test to assess the compatibility between
//	the IEI distribution associated to iei_sequence_A and iei_sequence_B.
//	The two sequences do not need to have the same number of elements.
//	The assessed p-value 'p' corresponds to the null-hypothesis of the two
//	sample distributions of IEI having the same parent distribution.
//	The K-S statistic 'd' is also provided.
//
//	This function is part of the SpiSeMe package.
//
//
//	REFERENCE:
//	The Kolmogorov-Smirnov test and its present implementation are thoroughly
//	described in W.H. Press, S.A. Teukolsky W.T. Vetterling and B.P. Flannery,
//	Numerical Recipes. The Art of Scientific Computing, 3rd Edition, 2007,
//	ISBN 0-521-88068-8.


#ifndef SPISEME_HEADER
	#define SPISEME_HEADER
	#include "spiSeMe.h"
#endif

double cdf_Q_KS(double);

spiSeMe_return_code spiSeMe_distribution_test(double & p, double & d, const std::vector<double> & iei_sequence_A, const std::vector<double> & iei_sequence_B)
{
	// --- Input parsing & validation
	unsigned int n_A = iei_sequence_A.size();
	unsigned int n_B = iei_sequence_B.size();
	if ((n_A < 2) || (n_B < 2)) {
		std::cerr << "ERROR (in spiSeMe_event_cross_correlation): input sequences are not actual vectors (length < 2).\n";
		return SSM_BAD_ARGUMENT;
	}
	for (int i = 0; i < n_A; i++) {
		if (iei_sequence_A[i] < 0) {
			std::cerr << "ERROR (in spiSeMe_distribution_test): invalid input sequence. One or more IEIs are negative.\n";
			return SSM_BAD_ARGUMENT;
		}
	}
	for (int i = 0; i < n_B; i++) {
		if (iei_sequence_B[i] < 0) {
			std::cerr << "ERROR (in spiSeMe_distribution_test): invalid input sequence. One or more IEIs are negative.\n";
			return SSM_BAD_ARGUMENT;
		}
	}

	// --- Key parameters
	double n_e = sqrt((double) n_A * (double) n_B / (double) (n_A + n_B));

	// --- K-S statistic: implementation following Numerical Recipes
	double d_statistic = 0.0;
	int j1 = 0, j2 = 0;
	double d1, d2, fn1 = 0.0, fn2 = 0.0;
	std::vector<double> data_A(iei_sequence_A);
	std::vector<double> data_B(iei_sequence_B);
	sort(data_A.begin(), data_A.end());
	sort(data_B.begin(), data_B.end());
	while (j1 < n_A && j2 < n_B) {
		d1 = data_A[j1];
		d2 = data_B[j2];
		if (d1 <= d2) {
			do
				fn1 = (double) ++j1 / (double) n_A;
			while (j1 < n_A && d1 == data_A[j1]);
		}
		if (d2 <= d1) {
			do
				fn2 = (double) ++j2 / (double) n_B;
			while (j2 < n_B && d2 == data_B[j2]);
		}
		if (fabs(fn2-fn1) > d_statistic)
			d_statistic = fabs(fn2-fn1);
	}
	d = d_statistic;
	double KS_statistic = (n_e + 0.12 + 0.11/n_e) * d_statistic;	// See Numerical Recipes Sec. 14.3.3

	p = cdf_Q_KS(KS_statistic);

	return SSM_SUCCESS;
}


double cdf_Q_KS(double z)
{
	double q, y;		// K-S distribution evaluated according to Sec. 6.14 of Numerical Recipes
	if (z == 0.) {
		q = 1.0;
	} else if (z < 1.18) {
		y = exp(-1.23370055013616983 / (z * z));
		q = 1.0 - 2.25675833419102515 * sqrt(-log(y)) * (y + pow(y, 9) + pow(y, 25) + pow(y, 49));
	} else {
		y = exp(-2.0 * z * z);
		q = 2.0 * (y - pow(y, 4) + pow(y, 9));
	}

	return q;
}
