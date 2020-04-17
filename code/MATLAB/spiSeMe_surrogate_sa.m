function ieiSurrogates = spiSeMe_surrogate_sa(ieiSequence, autocorrBinWidth, autocorrMaxLag, varargin)
% SPISEME_SURROGATE_SA Generates surrogates of Inter-Event-Intervals (IEI) sequences.
%
%	ieiSurrogates = SPISEME_SURROGATE_SA(ieiSequence, autocorrBinWidth, autocorrMaxLag, ...)
%	Generates, by means of the Simulated Annealing (SA) algorithm,
%	surrogate IEI sequences corresponding to the original
%	sequence ieiSequence.
%
%	The evaluation of the cost function relies on the assessment of
%	event autocorrelation. This assessment is implemented by the
%	function spiSeMe_event_autocorrelation within this package. The mandatory
%	parameters 'autocorrBinWidth', 'autocorrMaxLag' correspond to
%	the 'binWidth', 'maxLag' parameters of the spiSeMe_event_autocorrelation
%	function: please refer to its documentation for further details.
%
%	Options, passed as ('name',value) pairs:
%
%	'a'		Cooling factor (in references, alpha). Default 0.9.
%
%	'T'		Starting temperature. Default 0.1.
%
%	'C'		Target cost below which the routine is stopped. If
%			not specified, the function attempts to estimate a
%			reasonable target cost as the standard deviation of the
%			original sequence's autocorrelation. This estimate is not
%			necessarily a good choice for every sequence. Similarly,
%			different metrics (see 'costFunction') will lead to a
%			different accuracy under the default target cost.
%
%	'costFunction'	Sets the metric used to evaluate the cost function:
%			'max': Cost equals the maximim absolute difference
%			between original and surrogate autocorrelations
%			(this is the default setting).
%			'L1': Cost equals the sum over all bins of the absolute
%			differences between original and surrogate
%			autocorrelations.
%			'L2': Cost equals the sum over all bins of the
%			squared differences between original and surrogate
%			autocorrelations.
%
%	'nTotal'	Number of swaps to be carried out before cooling.
%			Default is N, where N is the sequence length.
%
%	'nSuccessful'	Number of accepted swaps to be carried out before
%			cooling. Default is N/2, where N is the sequence length.
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
%	The simulated annealing algorithm for surrogate generation was originally
%	proposed by T. Schreiber in Phys. Rev. Lett. 80 (1998), 2105,
%	<a href="matlab:web('https://doi.org/10.1103/PhysRevLett.80.2105')">doi:10.1103/PhysRevLett.80.2105</a>
%
%	Details on the autocorrelation assessment, as well as a comparative study
%	of surrogate generation algorithms, can be found in
%	Chaos 29 (2019), 121102, <a href="matlab:web('https://doi.org/10.1063/1.5138250')">doi:10.1063/1.5138250</a>
%
%	Citing this source is highly appreciated. Thank you.
%
%	See also spiSeMe_event_autocorrelation


	% --- Input parsing & validation
	ip = inputParser;
	addRequired(ip, 'ieiSequence', @isvector);
	addRequired(ip, 'autocorrBinWidth', @isnumeric);
	addRequired(ip, 'autocorrMaxLag', @isnumeric);
	addParameter(ip, 'a', 0.9, @isnumeric);
	addParameter(ip, 'T', 1.0, @isnumeric);
	addParameter(ip, 'C', -1, @isnumeric);
	validCosts = {'max', 'L1', 'L2'};
	checkCost = @(x) any(validatestring(x, validCosts));
	addParameter(ip, 'costFunction', 'max', checkCost);
	addParameter(ip, 'nTotal', -1, @isnumeric);
	addParameter(ip, 'nSuccessful', -1, @isnumeric);
	addParameter(ip, 'M', 1, @isnumeric);
	addParameter(ip, 'verbose', true, @islogical);
	ip.KeepUnmatched = false;
	parse(ip, ieiSequence, autocorrBinWidth, autocorrMaxLag, varargin{:});
	coolingFactor = ip.Results.a;
	startingTemperature = ip.Results.T;
	targetCost = ip.Results.C;
	if (~iscolumn(ieiSequence))
		ieiSequence = ieiSequence';
	end
	N = length(ieiSequence);
	M = fix(ip.Results.M);
	if (M < 1)
		error('Argument "M" must be a positive integer.')
	end
	beVerbose = ip.Results.verbose;
	if (autocorrBinWidth <= 0)
		error('Argument "autocorrBinWidth" must be a positive scalar value.');
	end
	if (autocorrMaxLag <= 0)
		error('Argument "autocorrMaxLag" must be a positive scalar value.');
	end
	if (autocorrMaxLag <= autocorrBinWidth)
		error('Invalid parameter: bin width must be smaller than maximum lag.')
	end
	if (startingTemperature < 0)
		error('Invalid parameter: starting temperature must be positive.')
	end
	if ((coolingFactor >= 1.0) || (coolingFactor <= 0))
		print('Invalid parameter: cooling factor must be positive and less than unity.')
	end
	if (strcmp(ip.Results.costFunction, 'max'))
		costIdx = 0;
	elseif (strcmp(ip.Results.costFunction, 'L1'))
		costIdx = 1;
	elseif (strcmp(ip.Results.costFunction, 'L2'))
		costIdx = 2;
	else
		error('Function argument "costFunction" does not match any of the allowed values ("max", "L1", "L2").');
	end

	% --- Assess defaults for omitted parameters
	nTotal = fix(ip.Results.nTotal);
	nSuccessful = fix(ip.Results.nSuccessful);
	if (nTotal <= 0)
		nTotal = N;
	end
	if (nSuccessful <= 0)
		nSuccessful = fix(N / 2);
	end

	% --- Assess autocorrelation of original sequence
	[originalAutocorr, ~] = spiSeMe_event_autocorrelation(ieiSequence, autocorrBinWidth, autocorrMaxLag);

	% --- If not assigned, estimate the target cost as the standard deviation of the Autocorrelation
	if (targetCost <= 0)
		targetCost = std(originalAutocorr);
		warning('Automatic setting of target cost might be too small for the routine to converge in a reasonable time.');
	end

	% --- Provide some information
	if (beVerbose)
		fprintf('\n### Starting simulated annealing routine ###\n');
		fprintf('###\tTarget cost: %f\n', targetCost);
		fprintf('###\tCooling factor: %f\n', coolingFactor);
		fprintf('###\tStarting T: %f\n', startingTemperature);
		fprintf('###\tSuccessful swaps before cooling: %d\n', nSuccessful);
		fprintf('###\tTotal swaps before cooling: %d\n\n', nTotal);
	end

	% --- Initialize output array (each of the M generated surrogates is a column)
	ieiSurrogates = [];

	% --- Generate surrogates: iterate M times
	for iter=1:M

		if(beVerbose)
			fprintf('# Surrogate number %d out of %d.\n\n', iter, M);
		end

		% --- Starting conditions
		rng('shuffle', 'twister');
		T = startingTemperature;
		nSwaps = 0;
		runCooling = true;
		runRandomization = true;
		sequenceSurrogateIEI(:,1) = ieiSequence(randperm(N));
		[surrogateAutocorr, ~] = spiSeMe_event_autocorrelation(sequenceSurrogateIEI, autocorrBinWidth, autocorrMaxLag);
		previousCost = computeCostFunction();

		% --- Outermost iteration: cooling stages
		while (runCooling)

			Ns = 0;
			Nt = 0;

			% --- Innermost iteration: swaps & Metropolis step
			while (runRandomization)

				indexesToBeSwapped = randi(N, 2, 1);
				sequenceSurrogateIEI(flip(indexesToBeSwapped)) = sequenceSurrogateIEI(indexesToBeSwapped);

				[surrogateAutocorr, ~] = spiSeMe_event_autocorrelation(sequenceSurrogateIEI, autocorrBinWidth, autocorrMaxLag);
				permutationCost = computeCostFunction();

				% --- Metropolis step
				deltaE = permutationCost - previousCost;
				if ((deltaE > 0) && (rand(1) > exp(-deltaE/T)))
					sequenceSurrogateIEI(flip(indexesToBeSwapped)) = sequenceSurrogateIEI(indexesToBeSwapped);
				else
					previousCost = permutationCost;
					Ns = Ns + 1;
				end
				Nt = Nt + 1;

				if ((previousCost <= targetCost) || (Nt == nTotal) || (Ns == nSuccessful))
					break;
				end
			end
			nSwaps = nSwaps + Nt;

			% --- Provide some information
			if (beVerbose)
				fprintf('Temperature: %f\tAccepted swaps: %d, Total swaps: %d\tCost: %f\n', T, Ns, Nt, previousCost);
			end

			% --- Cool by cooling factor
			T = coolingFactor * T;
			if (previousCost <= targetCost)
				break;
			end

		end

		if (beVerbose)
			fprintf('Target cost reached after %d total swaps.\n\n', nSwaps);
		end

		% --- Append sequence to the set of those already generated
		ieiSurrogates = [ieiSurrogates, sequenceSurrogateIEI];

	end

	% --- Cost function definition
	function c = computeCostFunction()
		if (costIdx == 0)
			c = max(abs(originalAutocorr - surrogateAutocorr));
		elseif (costIdx == 1)
			c = sum(abs(originalAutocorr - surrogateAutocorr));
		elseif (costIdx == 2)
			c = sum(abs(originalAutocorr - surrogateAutocorr).^2);
		end
	end

end
