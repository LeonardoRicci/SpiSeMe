""" Import libraries and package """
import numpy as np
import matplotlib.pyplot as plt
import spiSeMe.spiSeMe_surrogate_jodi as sjodi
import spiSeMe.spiSeMe_event_cross_correlation as scc

""" Load data """
# --- Read data from an ASCII file
iei_original_sequence_A = np.loadtxt('../data/iei_henon.dat')
iei_original_sequence_B = np.loadtxt('../data/iei_henon_bis.dat')

""" Display an excerpt of the original sequence """
time_interval_to_display = 50.0
arrival_times_original_A = np.cumsum(iei_original_sequence_A)
arrival_times_original_A = arrival_times_original_A[np.argwhere(arrival_times_original_A < time_interval_to_display)]
arrival_times_original_B = np.cumsum(iei_original_sequence_B)
arrival_times_original_B = arrival_times_original_B[np.argwhere(arrival_times_original_B < time_interval_to_display)]
fig, axes = plt.subplots(2, 1)
axes[0].stem(arrival_times_original_A, np.ones_like(arrival_times_original_A), linefmt='C0', markerfmt='C0o', basefmt='C0-')
axes[0].set_xticks([])
axes[0].set_yticks([])
axes[1].stem(arrival_times_original_B, np.ones_like(arrival_times_original_B), linefmt='C0', markerfmt='C0o', basefmt='C0-')
axes[1].set_xticks([])
axes[1].set_yticks([])
plt.xticks(np.arange(0,time_interval_to_display+1,time_interval_to_display/5))
plt.xlabel('Time')
fig.tight_layout()
plt.show()

""" Compute cross correlation of the original sequences """
bin_width = 0.06;
max_lag = 10.0;
C, lags = scc.spiSeMe_event_cross_correlation(iei_original_sequence_A, iei_original_sequence_B, bin_width, max_lag)

""" Cross correlation significance estimation """
# --- Surrogate generation
M = 500
surrogate_sequences_A = sjodi.spiSeMe_surrogate_jodi(iei_original_sequence_A, M=M)
surrogate_sequences_B = sjodi.spiSeMe_surrogate_jodi(iei_original_sequence_B, M=M)
# --- Assessment of cross correlation for each pair of surrogates
for i in range(0, M):
	C_temp, lags = scc.spiSeMe_event_cross_correlation(surrogate_sequences_A[i], surrogate_sequences_B[i], bin_width, max_lag)
	if (i == 0):
		C_surrogates = C_temp
	else:
		C_surrogates = np.vstack((C_surrogates, C_temp))
# --- Assessment of significance level
# --- The 1% significance level corresponds to the 6th largest value of each bin
C_surrogates = np.sort(C_surrogates, axis=0)
level = np.copy(C_surrogates[-6])

""" Display cross correlation of original sequences and significance level """
plt.plot(lags, C, label='Cross correlation')
plt.fill_between(lags, 0, level, facecolor=(0.3, 0.3, 0.3), alpha=0.5, label='1% significance level')
plt.xlim((-10, 10))
plt.ylim((0, 300))
plt.legend()
plt.xlabel('Lag')
plt.ylabel('Bin counts')
plt.show()
