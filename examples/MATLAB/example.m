%% Clear workspace & add code to path
clear;

% --- Add to path (until session ends) the directory where "spiSeMe_*.m" are stored.
% --- If necessary, modify the path accordingly.
codePath = '../../code/MATLAB';
addpath(codePath);

%% Select the sequence in input
% --- Uncomment one among the following labels:
% seqLabel = 'poisson';
% seqLabel = 'heartbeatlike';
seqLabel = 'henon';

%% Load data
% --- If existing, use 'readmatrix' to read data from an ASCII file.
% --- In some (< R2019a) versions of matlab, 'readmatrix' is unavailable:
% --- dlmread is used instead.
b = exist('readmatrix');
if (b == 2)
	originalSequence = readmatrix(strcat('../data/iei_', seqLabel, '.dat'));
else
	originalSequence = dlmread(strcat('../data/iei_', seqLabel, '.dat'));
end

%% Display an excerpt of the original sequence
timeIntervalToDisplay = 50;
originalArrivalTimes = cumsum(originalSequence);
originalArrivalTimes = originalArrivalTimes(originalArrivalTimes < timeIntervalToDisplay);
figure('Name', 'Excerpt of original sequence');
stem([0; originalArrivalTimes], ones(length(originalArrivalTimes) + 1, 1), 'Marker', 'none');
yticks([]);
xlabel('Time');

%% Generate surrogate IEIs
% --- As an example, three surrogate sequences are generated
M = 3;

% --- Select the surrogate generation method
% --- JODI
surrogateSequences = spiSeMe_surrogate_jodi(originalSequence, 'M', M);

% --- IAAFT
surrogateSequences = spiSeMe_surrogate_iaaft(originalSequence, 'ExactlyPreserve', 'distribution', 'M', M);

% --- SA (requires some extra parameters)
if (strcmp(seqLabel, 'poisson'))
	binWidth = 0.03;
	maxLag = 10.5;
elseif (strcmp(seqLabel, 'heartbeatlike'))
	binWidth = 0.03;
	maxLag = 36.0;
elseif (strcmp(seqLabel, 'henon'))
	binWidth = 0.06;
	maxLag = 21.0;
end
targetCost = 0.64;
surrogateSequences = spiSeMe_surrogate_sa(originalSequence, binWidth, maxLag, 'T', 0.1, 'a', 0.9, 'C', targetCost, 'CostFunction', 'max', 'nTotal', 1000, 'nSuccessful', 100, 'M', M);

% --- DITHERING
surrogateSequences = spiSeMe_surrogate_dither(originalSequence, 'M', M, 'DitherDistribution', 'uniform');

%% Display an excerpt of the original and surrogate sequences
timeIntervalToDisplay = 50;
figure('Name', 'Excerpt of original and surrogate sequence');

originalArrivalTimes = cumsum(originalSequence);
originalArrivalTimes = originalArrivalTimes(originalArrivalTimes < timeIntervalToDisplay);
subplot(M+1, 1, 1);
stem([0; originalArrivalTimes], ones(length(originalArrivalTimes) + 1,1), 'Marker', 'none');
xticks([]);
yticks([]);

for i=1:M
	surrogateArrivalTimes = cumsum(surrogateSequences);
	surrogateArrivalTimes = surrogateArrivalTimes(surrogateArrivalTimes(:,i) < timeIntervalToDisplay, i);
	subplot(M+1, 1, i+1);
	h = stem([0; surrogateArrivalTimes], ones(length(surrogateArrivalTimes) + 1, 1), 'Marker', 'none');
	h.Color = [0.8500 0.3250 0.0980];
	xticks([]);
	yticks([]);
end
xlabel('Time');
xticks('auto');

%% Compute autocorrelation of original and surrogate sequence
% --- The following parameters are the same as those reported in
% --- Chaos 29, 121102 (2019); doi:10.1063/1.5138250
if (strcmp(seqLabel, 'poisson'))
	binWidth = 0.03;
	maxLag = 10.5;
elseif (strcmp(seqLabel, 'heartbeatlike'))
	binWidth = 0.03;
	maxLag = 36.0;
elseif (strcmp(seqLabel, 'henon'))
	binWidth = 0.06;
	maxLag = 21.0;
end

% --- Autocorrelation of original sequence
[originalAutocorrelation, lags] = spiSeMe_event_autocorrelation(originalSequence, binWidth, maxLag);
% --- As an example, the autocorrelation of the first surrogate sequence is computed.
[surrogateAutocorrelation, ~] = spiSeMe_event_autocorrelation(surrogateSequences(:,2), binWidth, maxLag);

%% Display autocorrelation of original and surrogate sequence
figure('Name', 'Autocorrelation of original and surrogate sequence');
plot(lags, originalAutocorrelation); hold on;
plot(lags, surrogateAutocorrelation, 'Color', [0.8500 0.3250 0.0980]);
legend('Original', 'Surrogate');
xlabel('Lag');
ylabel('Autocorrelation');

%% Save the surrogate sequences to a file
% --- If existing, use 'writematrix' to save data to an ASCII file (tab-delimited)
% --- In some (< R2019a) versions of matlab, 'writematrix' is unavailable:
% --- dlmwrite is used instead.
b = exist('writematrix');
if (b == 2)
	writematrix(surrogateSequences, strcat('surr_', seqLabel, '.dat'), 'Delimiter', 'tab');
else
	dlmwrite(strcat('surr_', seqLabel, '.dat'), surrogateSequences, '\t');
end
