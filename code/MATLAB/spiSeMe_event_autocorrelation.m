function [A, lags] = spiSeMe_event_autocorrelation(ieiSequence, binWidth, maxLag)
% SPISEME_EVENT_AUTOCORRELATION Computes the sample autocorrelation of Inter-Event-Intervals (IEI) sequences.
%
%	[A, lags] = SPISEME_EVENT_AUTOCORRELATION(ieiSequence, binWidth, maxLag)
%	Computes the sample autocorrelation of ieiSequence.
%	The autocorrelation assessment is carried out by binning the lag axis in
%	bins of width binWidth and up to a maximum lag value maxLag.
%	The returned array "A" contains the counts corresponding to each
%	bin, normalized by N*N*binWidth/T, where N is the number of IEIs
%	and T the sequence duration. The central lag value of each bin is
%	returned in the "lags" array.
%
%	Both the binWidth (bin width) and maxLag (maximum lag value for
%	which autocorrelation has to be assessed) have to be positive real
%	numbers. In order to provide a meaningful assessment, binWidth must
%	be smaller than maxLag.
%	If maxLag is not an integer multiple of binWidth, the actual
%	maximum lag considered in the assessment of autocorrelation is
%	binWidth*ceil(maxLag/binWidth).
%
%	This function is part of the SpiSeMe package.
%
%
%	REFERENCE:
%	Details on the autocorrelation assessment, including the rationale
%	behind the normalization coefficient, can be found in
%	Chaos 29 (2019), 121102, <a href="matlab:web('https://doi.org/10.1063/1.5138250')">doi:10.1063/1.5138250</a>
%
%	Citing this source is highly appreciated. Thank you.


	% --- Input parsing & validation
	ip = inputParser;
	addRequired(ip, 'ieiSequence', @isvector);
	addRequired(ip, 'binWidth', @isnumeric);
	addRequired(ip, 'maxLag', @isnumeric);
	ip.KeepUnmatched = false;
	parse(ip, ieiSequence, binWidth, maxLag);
	if (~iscolumn(ieiSequence))
		ieiSequence = ieiSequence';
	end
	N = length(ieiSequence);
	if (binWidth <= 0)
		error('Argument "binWidth" must be a positive scalar value.')
	end
	if (maxLag <= 0)
		error('Argument "maxLag" must be a positive scalar value.')
	end
	if (maxLag <= binWidth)
		error('Invalid parameter: bin width must be smaller than maximum lag.')
	end
	if (any(ieiSequence < 0))
		error('Invalid input sequence: one or more IEIs are negative.')
	end

	% --- Compute arrival times out of IEIs
	arrivalTimes = [0; cumsum(ieiSequence)];
	sequenceDuration = arrivalTimes(length(arrivalTimes));

	% --- Initialize output arrays
	A = zeros(ceil(maxLag / binWidth), 1);
	nBins = length(A);
	lags = binWidth .* ([0:1:nBins-1]' + 0.5);

	% --- Compute autocorrelation
	for i = 2:1:length(arrivalTimes)
		deltaArrivalTimes = arrivalTimes(i) - arrivalTimes(1:i-1);
		idx = ceil(deltaArrivalTimes(deltaArrivalTimes < nBins*binWidth) / binWidth);
		% The maximum lag considered is nBins*binWidth, which is >= maxLag (see initializations)
	 	A(idx) = A(idx) + 1;
	end

	% --- Normalize autocorrelation [see Chaos 29, 121102 (2019)]
	A = A ./ (N * N * binWidth / sequenceDuration);

end
