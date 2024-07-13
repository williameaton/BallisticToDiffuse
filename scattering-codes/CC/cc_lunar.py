from calc_CC import calculate_CC
import numpy as np
import obspy as ob
import obspy.signal.filter as osf
import matplotlib.pyplot as plt
from wetools import norm, obspy_gen_mpl

# ------------------------------------------------------- EDIT INPUTS --------------------------------------------------

# Path to lunar data
data_path = f"../../data/lunar/filtered_removed_response/"

# Window length
T    = 50
fmin = 0.01
fmax = 2.0

# sub-matrix frequency band analysed separately -- see Fig 9 of paper (grey box)
submat_f1 = 0.6
submat_f2 = 1.1

# Uncomment out event you want to use:
#event = '1972'
#stn = '12'
#start = 575
#end = 5375

event = '1971b'
stn = '12'
start = 1780
end = 6580

#event = '1971a'
#stn = '14'
#start = 485
#end = 5285

# ----------------------------------------------------------------------------------------------------------------------

# Read in data
data = ob.read(f"{data_path}/{event}_{stn}_LP_lpz.mseed")


# Convert to numpy arrays
tr = data[0]
tr.data = norm(tr.data)
time, y = obspy_gen_mpl(tr)

# Generate figure & plot original time series
fig, ax = plt.subplots(3,1, figsize=(12,7))
ax[0].plot(time, tr.data, 'k')

# Slice the time-series and re-plot
tr = tr.slice(starttime=tr.stats.starttime + start  ,
          endtime=tr.stats.starttime + end)
time, y = obspy_gen_mpl(tr)
ax[0].plot(time + start, tr.data, 'r')

# Compute CC on the sliced time series
freq_out, C, C_mean, no_windows = calculate_CC(tr.data, dt=tr.stats.delta, T=T)

# Cutoff frequency bands:
indl = np.where(freq_out >= fmin)[0][0]
indh = np.where(freq_out >  fmax)[0][0]
freq_out = freq_out[indl:indh]

# Compute mean value for band-limited (fmin to fmax) CC matrix
C = C[indl:indh, indl:indh]
cm = np.mean(C[np.triu_indices_from(C, k=1)])
print("cm: ", cm)

# Plot the matrix with colourbar
im = ax[1].imshow(C, extent=(freq_out[0], freq_out[-1], freq_out[0], freq_out[-1]),origin='lower', cmap='magma_r')
fig.colorbar(mappable=im, ax=ax[1], use_gridspec=True)

# Add the submatrix box:
bandwidth = freq_out[-1]-freq_out[0]
ax[1].axhline(submat_f1,   xmin= (submat_f1-freq_out[0])/bandwidth, xmax = (submat_f2-freq_out[0])/bandwidth)
ax[1].axhline(submat_f2,   xmin= (submat_f1-freq_out[0])/bandwidth, xmax = (submat_f2-freq_out[0])/bandwidth)
ax[1].axvline(submat_f1,   ymin= (submat_f1-freq_out[0])/bandwidth, ymax = (submat_f2-freq_out[0])/bandwidth)
ax[1].axvline(submat_f2,   ymin= (submat_f1-freq_out[0])/bandwidth, ymax = (submat_f2-freq_out[0])/bandwidth)
ax[0].set_title(f"Mean CC: {cm}")


# Compute mean CC in submatrix:
indl = np.where(freq_out>=submat_f1)[0][0]
indh = np.where(freq_out <= submat_f2)[0][-1]
freq_sub = freq_out[indl:indh]
C_submatrix = C[indl:indh+1, indl:indh+1]

# Plot submatrix on separate axis
im = ax[2].imshow(C_submatrix, extent=(freq_sub[0], freq_sub[-1], freq_sub[0], freq_sub[-1]),origin='lower', cmap='magma_r')
cm_submatrix = np.mean(C_submatrix[np.triu_indices_from(C_submatrix, k=1)])
print("CMEAN of submatrix: ", cm_submatrix)

#plt.show()
# Save figure
plt.savefig(fname=f'./figs/JREV_cc_{event}_{stn}_T{T}.pdf', dpi='figure', format='pdf')