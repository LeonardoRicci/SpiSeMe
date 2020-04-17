function ieiSurrogates = spiSeMe_surrogate_iaaft(ieiSequence, varargin)
% SPISEME_SURROGATE_IAAFT Generates surrogates of Inter-Event-Intervals (IEI) sequences.
%
%	ieiSurrogates =  SPISEME_SURROGATE_IAAFT(ieiSequence, ...)
%	Generates, by means of the Iterative Amplitude Adjusted Fourier Transform
%	(IAAFT) algorithm, surrogate IEI sequences corresponding to the original
%	sequence ieiSequence.
%
%	Options, passed as ('name',value) pairs:
%
%	'exactlyPreserve' Specifies which function has to be preserved
%			exactly by the surrogate sequence. Can either be
%			'distribution' (default) or 'spectrum'.
%
%	'M'		Number of surrogate sequences to generate.
%			By default, M = 1. Each of the M surrogate sequences
%			is a column in the returned array, which therefore
%			has a size N x M, where N is the original sequence
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
%	The IAAFT (Iterative Amplitude-Adjusted Fourier Transform) algorithm
%	for surrogate generation was originally proposed by T. Schreiber
%	and A. Schmitz in Phys. Rev. Lett. 77 (1996), 635,
%	<a href="matlab:web('https://doi.org/10.1103/PhysRevLett.77.635')">doi:10.1103/PhysRevLett.77.635</a>
%
%	A comparative study of surrogate generation algorithms can be found in
%	Chaos 29 (2019), 121102, <a href="matlab:web('https://doi.org/10.1063/1.5138250')">doi:10.1063/1.5138250</a>
%
%	Citing this source is highly appreciated. Thank you.


	% --- Input parsing & validation
	ip = inputParser;
	addRequired(ip, 'ieiSequence', @isvector);
	validExactlyPreserve = {'distribution', 'spectrum'};
	checkExactlyPreserve = @(x) any(validatestring(x, validExactlyPreserve));
	addParameter(ip, 'exactlyPreserve', 'distribution', checkExactlyPreserve);
	addParameter(ip, 'M', 1, @isnumeric);
	addParameter(ip, 'verbose', true, @islogical);
	ip.KeepUnmatched = false;
	parse(ip, ieiSequence, varargin{:});
	preserve = ip.Results.ExactlyPreserve;
	if (~iscolumn(ieiSequence))
		ieiSequence = ieiSequence';
	end
	N = length(ieiSequence);
	M = fix(ip.Results.M);
	if (M < 1)
		error('Function argument "M" must be a positive integer.')
	end
	beVerbose = ip.Results.verbose;

	% --- Provide some information
	if (beVerbose)
		fprintf('\n### Starting IAAFT routine ###\n');
		fprintf('###\tExactly preserving: %s\n', preserve);
	end

	% --- Assess spectrum & distribution of original sequence
	spectrumMagnitudeOriginal = abs(fft(ieiSequence));
	distributionOriginal = sort(ieiSequence);

	% --- Initialize output array (each of the M generated surrogates is a column)
	ieiSurrogates = [];

	% --- Generate surrogates: iterate M times
	for iter=1:M
		if(beVerbose)
			fprintf('# Surrogate number %d out of %d.\n', iter, M);
		end

		% --- Starting conditions
		rng('shuffle', 'twister');
		runIterations = true;
		sequenceSurrogateIEI(:,1) = ieiSequence(randperm(N));
		lastSequenceSurrogateIEI = sequenceSurrogateIEI;

		% --- Iterative algorithm (switched depending on target function to be matched exactly.
		nIter = 0;
		while (runIterations)

			nIter = nIter + 1;

			spectrumSurrogate = fft(sequenceSurrogateIEI);
			phasesSurrogate = angle(spectrumSurrogate);
			spectrumSurrogate = spectrumMagnitudeOriginal .* (cos(phasesSurrogate) + 1i.*sin(phasesSurrogate));
			sequenceSurrogateIEI = ifft(spectrumSurrogate, 'symmetric');

			[~, iSorted] = sort(sequenceSurrogateIEI);
			sequenceSurrogateIEI(iSorted) = distributionOriginal;

			if (sequenceConverged())
				if (strcmp(preserve, 'spectrum'))
					spectrumSurrogate = fft(sequenceSurrogateIEI);
					phasesSurrogate = angle(spectrumSurrogate);
					spectrumSurrogate = spectrumMagnitudeOriginal .* (cos(phasesSurrogate) + 1i.*sin(phasesSurrogate));
					sequenceSurrogateIEI = ifft(spectrumSurrogate, 'symmetric');
				end
				break;
			end

			lastSequenceSurrogateIEI = sequenceSurrogateIEI;
		end

		if (beVerbose)
			fprintf('Process converged after %d iterations.\n\n', nIter);
		end

		% --- Append sequence to the set of those already generated
		ieiSurrogates = [ieiSurrogates, sequenceSurrogateIEI];
	end

	% --- Define check of convergence (sequence did not change with respect to last iteration)
	function convergence = sequenceConverged()
		c = sum(abs(lastSequenceSurrogateIEI - sequenceSurrogateIEI));
		if (c == 0)
			convergence = true;
		else
			convergence = false;
		end
	end

end
