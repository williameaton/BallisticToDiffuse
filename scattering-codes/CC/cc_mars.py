from calc_CC import calculate_CC
import numpy as np
import matplotlib.pyplot as plt


# Window size:
T = 50

# Define frequency bandwidth
fmin = 0.1
fmax = 2.0

# Datapath, event and timestep:
data_path = f"../../data/mars/deglitched_downsampled/"
event='S0235'
dt = '0.151'

# Load data that is deglitched and downsampled from user-defined data_path
data = np.loadtxt(f"{data_path}/HREV_NEW_{event}_delta_{dt}_{fmin}_{fmax}")

# Get number of channels:
channels = np.shape(data)[0]

# Get time
time = np.linspace(0, len(data[0,:])*float(dt),  len(data[0,:]))


# Generate figure
fig, ax = plt.subplots(2,channels)

# Loop channels:
for chl in range(channels):
    # Plot original time series
    ax[0,chl].plot(time, data[chl,:])

    # Compute CC matrix
    freq_out, C, C_mean, no_windows = calculate_CC(data[chl,:], dt=float(dt), T=T)

    # Slice CC matrix based on frequencies
    indl = np.where(freq_out >= fmin)[0][0]
    indh = np.where(freq_out >  fmax)[0][0]
    freq_out = freq_out[indl:indh]
    C = C[indl:indh, indl:indh]

    # Compute mean CC within frequency band [fmin, fmax]
    cm = np.mean(C[np.triu_indices_from(C, k=1)])
    print("cm: ", cm)

    # Plot CC matrix and colourbar
    im = ax[1,chl].imshow(C, extent=(freq_out[0], freq_out[-1], freq_out[0], freq_out[-1]),origin='lower', cmap='magma_r')
    fig.colorbar(mappable=im, ax=ax[1,chl], use_gridspec=True)
    ax[0,chl].set_title(f"Mean CC: {cm}")


plt.show()

# Save figure
plt.savefig(fname=f'./figs/REVIEW_cc_{event}_{dt}_T{T}.pdf', dpi='figure', format='pdf')