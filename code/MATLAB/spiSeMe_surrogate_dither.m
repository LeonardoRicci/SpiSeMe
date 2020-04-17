function ieiSurrogates = spiSeMe_surrogate_dither(ieiSequence, varargin)
% SPISEME_SURROGATE_DITHER Generates surrogates of Inter-Event-Intervals (IEI) sequences.
%
%	ieiSurrogates = SPISEME_SURROGATE_DITHER(ieiSequence, ...)
%	Generates, by means of dithering, surrogate IEI sequences corresponding
%	to the original sequence ieiSequence.
%
%	Dithering consists in adding random time intervals to the time
%	coordinate of each event in the sequence. Random intervals follow a
%	zero-mean distribution whose shape and width can be set by the
%	user.
%
%	Options, passed as ('name',value) pairs:
%
%	'ditherDistribution' Specifies from which distribution have the dithering
%			elements to be drawn. Can either be 'uniform'
%			(default), 'normal', or 'triangular'.
%
%	'distributionParameter'	(Approximate) half-width of the distribution from
%			which dithering elements have to be drawn (see
%			'Distribution'). In the case of 'uniform', the
%			distribution is given by U(-D, D). In the
%			case of 'normal', the distribution has standard
%			deviation equal to D/5.0. In the case of 'triangular',
%			the distribution is the symmetric triangular distribution
%			T(-D, D). By default, D is estimated as half the minimum
%			IEI within the original sequence.
%
%	'M'		Number of surrogate sequences to be generated. By
%			default, M = 1. Each of the M surrogate sequences
%			is a column in the returned array, which therefore
%			has size N x M, where N is the original sequence
%			length.
%
%	'verbose'	Sets the verbosity of the function. If true (default),
%			all messages are printed on the command line. If
%			false, only critical errors are displayed.
%
%	This function is part of the SpiSeMe package.
%
%
%	REFERENCE:
%	A comparative study of surrogate generation algorithms can be found in
%	Chaos 29 (2019), 121102, <a href="matlab:web('https://doi.org/10.1063/1.5138250')">doi:10.1063/1.5138250</a>
%
%	Citing this source is highly appreciated. Thank you.


	% --- Input parsing & validation
	ip = inputParser;
	addRequired(ip, 'ieiSequence', @isvector);
	validDistributions = {'uniform', 'normal', 'triangular'};
	checkDistribution = @(x) any(validatestring(x, validDistributions));
	addParameter(ip, 'ditherDistribution', 'uniform', checkDistribution);
	addParameter(ip, 'distributionParameter', 0.0, @isnumeric);
	addParameter(ip, 'M', 1, @isnumeric);
	addParameter(ip, 'verbose', true, @islogical);
	ip.KeepUnmatched = false;
	parse(ip, ieiSequence, varargin{:});
	if (~iscolumn(ieiSequence))
		ieiSequence = ieiSequence';
	end
	M = fix(ip.Results.M);
	if (M < 1)
		error('Function argument "M" must be a positive integer.')
	end
	beVerbose = ip.Results.verbose;
	ditherDistribution = ip.Results.DitherDistribution;
	DistributionParameter = ip.Results.DistributionParameter;
	if ((~strcmp(ditherDistribution, 'uniform')) && (~strcmp(ditherDistribution, 'normal')) && (~strcmp(ditherDistribution, 'triangular')))
		error('Function argument "distribution" does not match any of the allowed values ("uniform", "normal", "triangular").');
	end

	% --- If not provided, estimate DistributionParameter parameter
	if (DistributionParameter <= 0)
		DistributionParameter = min(ieiSequence)/2.0;
	end

	% --- Provide some information
	if (beVerbose)
		fprintf('\n### Starting dithering routine ###\n');
		fprintf('###\tDither distribution: %s\n', ditherDistribution);
		if (strcmp(ditherDistribution, 'uniform'))
			fprintf('###\t\tU(-D,D), with D = %f\n', DistributionParameter);
		elseif (strcmp(ditherDistribution, 'normal'))
			fprintf('###\t\tN(0,s^2), with s = %f\n', DistributionParameter/5.0);
		elseif (strcmp(ditherDistribution, 'triangular'))
			fprintf('###\t\tT(-D,D), with D = %f\n', DistributionParameter);
		end
	end

	% --- Build arrival times to be dithered
	originalArrivalTimes = [0; cumsum(ieiSequence)];

	% --- Initialize output array (each of the M generated surrogates is a column)
	ieiSurrogates = [];

	% --- Generate surrogates: iterate M times
	for iter=1:M

		if (beVerbose)
			fprintf('# Surrogate number %d out of %d.\n', iter, M);
		end

		% --- Initialize surrogate ranks generation
		rng('shuffle', 'twister');
		arrivalTimes = originalArrivalTimes;

		% --- Apply dithering to arrival times
		if (strcmp(ditherDistribution, 'uniform'))
			deltas = DistributionParameter .* (-1.0 + 2.0 .* rand(length(arrivalTimes), 1));
			arrivalTimes = arrivalTimes + deltas;
		elseif (strcmp(ditherDistribution, 'normal'))
			deltas = (DistributionParameter / 5.0) .* randn(length(arrivalTimes), 1);
			arrivalTimes = arrivalTimes + deltas;
		elseif (strcmp(ditherDistribution, 'triangular'))
			deltas = DistributionParameter .* (rand(length(arrivalTimes), 1) - rand(length(arrivalTimes), 1));
			arrivalTimes = arrivalTimes + deltas;
		end

		% --- Transform arrival times back into IEIs
		thisSurrogate = arrivalTimes(2:length(arrivalTimes)) - arrivalTimes(1:length(arrivalTimes)-1);

		if (beVerbose)
			fprintf('Process ended.\n\n');
		end

		% --- Append sequence to the set of those already generated
		ieiSurrogates = [ieiSurrogates, thisSurrogate];
	end
end
