%% Clear workspace & add code to path
clear;

% --- Add to path (until session ends) the directory where "spiSeMe_*.m" are stored.
% --- If necessary, modify the path accordingly.
codePath = '../../code/MATLAB';
addpath(codePath);

%% Load data
% --- If existing, use 'readmatrix' to read data from an ASCII file.
% --- In some (< R2019a) versions of matlab, 'readmatrix' is unavailable:
% --- dlmread is used instead.
b = exist('readmatrix');
if (b == 2)
	sequenceA = readmatrix(strcat('../data/iei_henon.dat'));
	sequenceB = readmatrix(strcat('../data/iei_henon_bis.dat'));
else
	sequenceA = dlmread(strcat('../data/iei_henon.dat'));
	sequenceB = dlmread(strcat('../data/iei_henon_bis.dat'));
end

%% Display an excerpt of the original sequences
timeIntervalToDisplay = 50;
arrivalTimesA = cumsum(sequenceA);
arrivalTimesA = arrivalTimesA(arrivalTimesA < timeIntervalToDisplay);
arrivalTimesB = cumsum(sequenceB);
arrivalTimesB = arrivalTimesB(arrivalTimesB < timeIntervalToDisplay);
figure('Name', 'Excerpt of original sequences');
subplot(2, 1, 1);
stem([0; arrivalTimesA], ones(length(arrivalTimesA) + 1, 1), 'Marker', 'none');
yticks([]);
subplot(2, 1, 2);
stem([0; arrivalTimesB], ones(length(arrivalTimesB) + 1, 1), 'Marker', 'none');
yticks([]);
xlabel('Time');

%% Compute cross correlation of the original sequences
binWidth = 0.06;
maxLag = 10.0;
[C, lags] = spiSeMe_event_cross_correlation(sequenceA, sequenceB, binWidth, maxLag);

%% Cross correlation significance estimation
% --- Surrogate generation
M = 500;
surrogateSequencesA = spiSeMe_surrogate_jodi(sequenceA, 'M', M);
surrogateSequencesB = spiSeMe_surrogate_jodi(sequenceB, 'M', M);
% --- Assessment of cross correlation for each pair of surrogates
Csurr = [];
for i = 1:M
	[Ctemp, ~] = spiSeMe_event_cross_correlation(surrogateSequencesA(:,i), surrogateSequencesB(:,i), binWidth, maxLag);
	Csurr = [Csurr, Ctemp];
end
% --- Assessment of significance level
% --- The 1% significance level corresponds to the 6th largest value of each bin
Csurr = sort(Csurr, 2, 'descend');
level = Csurr(:,6);

%% Display cross correlation of original sequences and significance level
figure('Name', 'Cross correlation assessment');
plot(lags, C); hold on;
area(lags, level, 'FaceColor', [0.3 0.3 0.3], 'FaceAlpha', 0.5, 'LineStyle', 'none');
legend('Cross correlation', '1% significance level');
xlim([-10 10])
ylim([0 300])
xlabel('Lag');
ylabel('Bin counts');
