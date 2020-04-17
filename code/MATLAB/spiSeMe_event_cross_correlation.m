function [C, lags] = spiSeMe_event_cross_correlation(ieiSequenceA, ieiSequenceB, binWidth, maxLag)
% SPISEME_EVENT_CROSS_CORRELATION Computes the sample cross-correlation between
%	two Inter-Event-Intervals (IEI) sequences.
%
%	[C, lags] = SPISEME_EVENT_CROSS_CORRELATION(ieiSequenceA, ieiSequenceB, binWidth, maxLag)
%	Computes the sample cross correlation between ieiSequenceA and ieiSequenceB.
%	The cross correlation assessment is carried out by binning the lag axis in
%	bins of width binWidth for lag values from -max_lag to max_lag.
%	The returned array "C" contains the counts corresponding to each
%	bin. The central lag value of each bin is returned in the "lags" array.
%
%	Both the binWidth (bin width) and maxLag (maximum lag value for
%	which autocorrelation has to be assessed) have to be positive real
%	numbers. In order to provide a meaningful assessment, binWidth must
%	be smaller than maxLag.
%	If maxLag is not an half-integer multiple of binWidth, the actual
%	maximum lag considered in the assessment of autocorrelation is
%	binWidth*(ceil(maxLag/binWidth) + 0.5).
%
%	This function is part of the SpiSeMe package.


	% --- Input parsing & validation
	ip = inputParser;
	addRequired(ip, 'ieiSequenceA', @isvector);
	addRequired(ip, 'ieiSequenceB', @isvector);
	addRequired(ip, 'binWidth', @isnumeric);
	addRequired(ip, 'maxLag', @isnumeric);
	ip.KeepUnmatched = false;
	parse(ip, ieiSequenceA, ieiSequenceB, binWidth, maxLag);
	if (~iscolumn(ieiSequenceA))
		ieiSequenceA = ieiSequenceA';
	end
	if (~iscolumn(ieiSequenceB))
		ieiSequenceB = ieiSequenceB';
	end
	if (binWidth <= 0)
		error('Argument "binWidth" must be a positive scalar value.')
	end
	if (maxLag <= 0)
		error('Argument "maxLag" must be a positive scalar value.')
	end
	if (maxLag <= binWidth)
		error('Invalid parameter: bin width must be smaller than maximum lag.')
	end
	if (any(ieiSequenceA < 0) || any(ieiSequenceB < 0))
		error('Invalid input sequence: one or more IEIs are negative.')
	end

	% --- Compute arrival times out of IEIs
	arrivalTimesA = [0; cumsum(ieiSequenceA)];
	arrivalTimesB = [0; cumsum(ieiSequenceB)];

	% --- Initialize output arrays
	nBins = 2 * ceil(maxLag / binWidth) + 1;
	C = zeros(nBins, 1);
	lags = binWidth .* [-floor(nBins/2):1:floor(nBins/2)]';
	actualMaxLag = abs(lags(1)) + 0.5*binWidth;

	% --- Compute cross correlation
	for i = 1:1:length(arrivalTimesA)
		deltaArrivalTimes = arrivalTimesA(i) - arrivalTimesB(:);
		idx = floor(deltaArrivalTimes(abs(deltaArrivalTimes) < actualMaxLag)/binWidth + 0.5) + floor(nBins/2) + 1;
		for j = 1:length(idx)
			C(idx(j)) = C(idx(j)) + 1;
		end
	end

end
