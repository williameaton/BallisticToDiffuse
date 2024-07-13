# Generate example 'made up' seismogram used for Apx C (Apx Fig 2)
import matplotlib.pyplot as plt
import numpy as np
import colorednoise as cn
from wetools import norm
from obspy.core.trace import Trace

beta = 0                    # white noise
samples = int(1000 / 0.15)  # number of samples to generate - means we have dt of 0.005


# Create original sinusoid
x = np.linspace(0, 1000, samples)
sinusoid = norm(0.5 * np.sin(2 * np.pi * x) + np.sin(3 * np.pi * x) + 0.7 * np.sin(3.6 * np.pi * x))

# Set number of coefficients between 0 and 300
num_coeffs = 100
max_C = 300
C = np.linspace(0, max_C, num_coeffs)


vals = [[], [], [], []]

# Generate noise
n1 = norm(cn.powerlaw_psd_gaussian(beta, samples))
N = np.random.normal(loc=0, scale=1, size=samples)

clrs = ['r', 'b', 'g', 'orange']

# Coefficient of different effects (see apx text)
precursor_coeff     = 250
decay_coeff         = 250
taper_noise_coeff   = 250
white_noise_coeff   = 250

# Precursor effect
front_width = precursor_coeff / 15
taper = np.zeros(samples)
taper[np.logical_and(x >= 5, x < 5 + front_width)] = (1 - np.cos(
    (x[np.logical_and(x >= 5, x < 5 + front_width)] - 5) * np.pi / (front_width))) / 2
taper[np.logical_and(x >= 5 + front_width, x < 6 + 2 * front_width)] = 1


# Coda Q effect:
taper[x >= 6 + 2 * front_width] = norm(np.exp(-x[x >= 6 + 2 * front_width] / (1 + decay_coeff)))
# taper[x >= 6+ 2*front_width] = norm(np.exp(-x[x >= 6+ 2*front_width] /(1+ 50))  )


# Noise in taper
n2 = norm(np.power(np.abs(N), 5 - 5 * taper_noise_coeff / max_C))
taper[x >= 5] = norm(norm(taper[x >= 5] * n2[x >= 5]) * taper[x >= 5])


# Noise in trace:
timeseries = norm((sinusoid + 2 * n1 * white_noise_coeff/ max_C) * taper)
tr = Trace()
tr.data = timeseries
tr.stats.delta = x[2] - x[1]


# Filter timeseries
tr = tr.filter(type='bandpass', freqmin=0.01, freqmax=2.0)
timeseries = norm(tr.data)


# PLOT FIGURE:
fig, ax = plt.subplots(3, 2, figsize=(13, 7.5),sharex=True)
# Original sinusoid
ax[0,0].plot(x,norm(sinusoid), 'k')
ax[0,0].set_xlim([0,1000])
ax[0,0].set_ylim([-1,1])
ax[0,0].set_title('Original sinusoid')

# Taper with just decay
ax[0,1].plot(x, taper, 'k')

# Noise generated:
#ax[0,1].plot(x, n2,  clrs[i])
ax[0,1].set_xlim([0, 300])
ax[0,1].set_ylim([0, 1])
ax[0,1].set_title('Generated white noise')

# Actual generated seismogram:
ax[1,0].plot(x, taper, 'k', alpha=0.4)
ax[1,1].plot(x, timeseries,  'k')
ax[1,0].set_xlim([0, 300])
ax[1,0].set_title('Resulting seismogram')

# Save output:
plt.savefig(fname=f'example_seis_mseapx.pdf', dpi='figure', format='pdf')


