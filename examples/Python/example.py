""" Import libraries and package """
import numpy as np
import matplotlib.pyplot as plt
import spiSeMe.spiSeMe_surrogate_jodi as sjodi
import spiSeMe.spiSeMe_surrogate_iaaft as siaaft
import spiSeMe.spiSeMe_surrogate_sa as ssa
import spiSeMe.spiSeMe_surrogate_dither as sdit
import spiSeMe.spiSeMe_event_autocorrelation as sac


""" Select the sequence in input """
# --- Uncomment one among the following labels:
#seq_label = 'poisson'
#seq_label = 'heartbeatlike'
seq_label = 'henon'

""" Load data """
# --- Read data from an ASCII file
iei_original_sequence = np.loadtxt('../data/iei_' + seq_label + '.dat')

""" Display an excerpt of the original sequence """
time_interval_to_display = 50.0
arrival_times_original = np.cumsum(iei_original_sequence)
arrival_times_original = arrival_times_original[np.argwhere(arrival_times_original < time_interval_to_display)]
plt.figure('Excerpt of original sequence')
plt.stem(arrival_times_original, np.ones_like(arrival_times_original), linefmt='C0', markerfmt='C0o', basefmt='C0-')
plt.xlabel('Time')
plt.yticks([])
plt.show()
plt.close()

""" Generate surrogate IEIs """
# --- As an example, three surrogate sequences are generated
M = 3

# --- JODI
iei_surrogate_sequences = sjodi.spiSeMe_surrogate_jodi(iei_original_sequence, M)

# --- IAAFT
iei_surrogate_sequences = siaaft.spiSeMe_surrogate_iaaft(iei_original_sequence, exactly_preserve='distribution', M=M);

# --- SA (requires some extra parameters)
if (seq_label == 'poisson'):
	bin_width = 0.03;
	max_lag = 10.5;
elif (seq_label == 'heartbeatlike'):
	bin_width = 0.03;
	max_lag = 36.0;
elif (seq_label == 'henon'):
	bin_width = 0.06;
	max_lag = 21.0;
target_cost = 0.64
iei_surrogate_sequences = ssa.spiSeMe_surrogate_sa(iei_original_sequence, bin_width, max_lag, T=0.1, a=0.9, C=target_cost, cost_function='max', n_total=1000, n_successful=100, M=M);

# --- DITHERING
iei_surrogate_sequences = sdit.spiSeMe_surrogate_dither(iei_original_sequence, dither_distribution='uniform', M=M);

""" Display an excerpt of the original and surrogate sequences """
time_interval_to_display = 50.0
fig, axes = plt.subplots(M + 1, 1)
fig.tight_layout()
axes[0].stem(arrival_times_original, np.ones_like(arrival_times_original), linefmt='C0', markerfmt='C0o', basefmt='C0-')
axes[0].set_xlabel('Time')
axes[0].set_yticks([])

for i in range(0, M):
	arrival_times_surrogate = np.cumsum(iei_surrogate_sequences[i])
	arrival_times_surrogate = arrival_times_surrogate[np.argwhere(arrival_times_surrogate < time_interval_to_display)]
	axes[i+1].stem(arrival_times_surrogate, np.ones_like(arrival_times_surrogate), linefmt='C1', markerfmt='C1o', basefmt='C1-')
	axes[i+1].set_xticks([])
	axes[i+1].set_yticks([])

plt.xticks(np.arange(0,51,10))
plt.xlabel('Time')
plt.show()

""" Compute autocorrelation of original and surrogate sequence """
# --- The following parameters are the same as those reported in
# --- Chaos 29, 121102 (2019); doi:10.1063/1.5138250
if (seq_label == 'poisson'):
	bin_width = 0.03
	max_lag = 10.5
elif (seq_label == 'heartbeatlike'):
	bin_width = 0.03
	max_lag = 36.0
elif (seq_label == 'henon'):
	bin_width = 0.06
	max_lag = 21.0
# --- Autocorrelation of original sequence
A_original, lags = sac.spiSeMe_event_autocorrelation(iei_original_sequence, bin_width, max_lag)
# --- As an example, the autocorrelation of the first surrogate sequence is computed.
A_surrogate, lags = sac.spiSeMe_event_autocorrelation(iei_surrogate_sequences[0], bin_width, max_lag)

""" Display autocorrelation of original and surrogate sequence """
plt.plot(lags, A_original, label='Original')
plt.plot(lags, A_surrogate, label='Surrogate')
plt.legend()
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.show()

# """ Save the surrogate sequence to a file """
numpy.savetxt('surr_' + seq_label + '.dat', np.transpose(iei_surrogate_sequences), fmt='%.4e', delimiter='\t')
