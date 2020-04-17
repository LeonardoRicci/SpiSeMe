function [p, d] = spiSeMe_distribution_test(ieiSequenceA, ieiSequenceB)
% SPISEME_DISTRIBUTION_TEST Tests compatibility of the IEI distribution of two sequences.
%
%	[p, d] = SPISEME_DISTRIBUTION_TEST(ieiSequenceA, ieiSequenceB)
%	Carries out a Kolmogorov-Smirnov test to assess the compatibility between
%	the IEI distribution associated to ieiSequenceA and ieiSequenceB.
%	The two sequences do not need to have the same number of elements.
%	The returned p-value 'p' corresponds to the null-hypothesis of the two
%	sample distributions of IEI having the same parent distribution.
%	The K-S statistic 'd' is also returned.
%
%	This function is part of the SpiSeMe package.
%
%
%	REFERENCE:
%	The Kolmogorov-Smirnov test and its present implementation are thoroughly
%	described in W.H. Press, S.A. Teukolsky W.T. Vetterling and B.P. Flannery,
%	Numerical Recipes. The Art of Scientific Computing, 3rd Edition, 2007,
%	ISBN 0-521-88068-8.


	% --- Input parsing & validation
	ip = inputParser;
	addRequired(ip, 'ieiSequenceA', @isvector);
	addRequired(ip, 'ieiSequenceB', @isvector);
	ip.KeepUnmatched = false;
	parse(ip, ieiSequenceA, ieiSequenceB);
	if (~iscolumn(ieiSequenceA))
		ieiSequenceA = ieiSequenceA';
	end
	if (~iscolumn(ieiSequenceB))
		ieiSequenceB = ieiSequenceB';
	end
	if (any(ieiSequenceA < 0))
		error('Invalid input sequence: one or more IEIs are negative.')
	end
	if (any(ieiSequenceB < 0))
		error('Invalid input sequence: one or more IEIs are negative.')
	end

	% --- Key parameters
	nA = length(ieiSequenceA);
	nB = length(ieiSequenceB);
	nE = sqrt(nA * nB / (nA + nB));

	% --- K-S statistic: implementation following Numerical Recipes
	d = 0.0;
	j1 = 1;
	j2 = 1;
	fn1 = 0.0;
	fn2 = 0.0;
	data_A = sort(ieiSequenceA);
	data_B = sort(ieiSequenceB);
	while (j1 <= nA && j2 <= nB)
		d1 = data_A(j1);
		d2 = data_B(j2);
		if (d1 <= d2)
			while (j1 <= nA && d1 == data_A(j1))
				j1 = j1 + 1;
				fn1 = j1 / nA;
			end
		end
		if (d2 <= d1)
			while (j2 <= nB && d2 == data_B(j2))
				j2 = j2 + 1;
				fn2 = j2 / nB;
			end
		end
		if (abs(fn2 - fn1) > d)
			d = abs(fn2 - fn1);
		end
	end
	KSstatistic = (nE + 0.12 + 0.11/nE) * d; % See Numerical Recipes Sec. 14.3.3

	p = qks(KSstatistic);

	% --- Definition of CDF to compute p-value of K-S statistics (see Numerical Recipes Sec. 6.14.12)
 	function q = qks(z)
		if (z == 0.)
			q = 1.0;
		elseif (z < 1.18)
			y = exp(-1.23370055013616983 / z^2);
			q = 1.0 - 2.25675833419102515 * sqrt(-log(y)) * (y + y^9 + y^25 + y^49);
		else
			y = exp(-2.0 * z^2);
			q = 2.0 * (y - y^4 + y^9);
		end
 	end

end
