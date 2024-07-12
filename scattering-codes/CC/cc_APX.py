import numpy as np
import matplotlib.pyplot as plt
import colorednoise as cn
from wetools import norm
from calc_CC import calculate_CC

# Create white noise
T = 10                                             # Time length of each window
dt   = 0.005                                       # Timestep
beta = 0                                           # Coloured noise beta coefficient                          # white noise
samples = int(250/dt)                              # Timeseries length for 250 second timeseries
n1 = norm(cn.powerlaw_psd_gaussian(beta, samples)) # Coloured noise generated
fmin = 0.01                                        # Minimum frequency
fmax = 2.0                                         # Maximum frequency



# Outputs of calculate CC:
#   frequency values
#   correlation matrix
#   mean of correlation matrix
#   number of windows used
freq_out, C, C_mean, no_windows = calculate_CC(n1, dt=dt, T=T)

# Slice out frequencies desired e.g. between 0.01 and 2 hz
indl = np.where(freq_out >= 0.01)[0][0]
indh = np.where(freq_out > 2)[0][0]
freq_out = freq_out[indl:indh]
C = C[indl:indh, indl:indh]

# Computes mean only within this sliced (frequency banded) matrix
cm = np.mean(C[np.triu_indices_from(C, k=1)])
print("cm: ", cm)

# Plot matrix:
# Since it is noise it should be ~ an identity matrix
fig, ax = plt.subplots()
im = ax.imshow(C, extent=(freq_out[0], freq_out[-1], freq_out[0], freq_out[-1]), origin='lower', cmap='magma_r')
fig.colorbar(mappable=im, ax=ax, use_gridspec=True)

plt.show()