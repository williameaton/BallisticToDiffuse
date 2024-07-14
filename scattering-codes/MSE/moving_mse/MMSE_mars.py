import matplotlib.gridspec as gs
import numpy as np
import pyentrp.entropy as ent
import matplotlib.pyplot as plt
from gen_MMSE_slices import gen_mmse_slices
from wetools import norm

# Input parameters:
event      = 'S0235'    # event name
dt         = '0.151'    # downsampled timestep for time series
channel    = 'Z'        # channel
slice_time = 50         # seconds
max_tau    = 100        # max scale
tolerance  = 0.03       # mse tolerance
vmax       = 3.5        # max value plot colourbar

# Scale array
m = np.arange(1, max_tau+1)

# Number of timesteps in single slice
slice_length = int(slice_time/float(dt))

# Load deglitched and downsampled data
data_path = f"../../../../data/mars/deglitched_downsampled/"
event_tr = np.loadtxt(f"{data_path}/FORMMSE_REV_NEW_{event}_delta_{dt}_0.1_2.0")

# Convert to np arrays:
trace = norm(event_tr[0,:])
length = len(trace)
time = np.linspace(0, float(dt)*length, length)


# Now need to generate slices for MMSE and time arrays:
timeslices, no_slices = gen_mmse_slices(trace        = time,
                                       slice_length = slice_length,
                                       overlap      = 0,
                                       normalise    = False)

slices, no_slices     = gen_mmse_slices(trace        = trace,
                                       slice_length = slice_length,
                                       overlap      = 0,
                                       normalise = True)

# ==================== FORMAT OUTPUT FIGURE: ======================
fig = plt.figure(constrained_layout=True, figsize=(13,6.5),)
spec = gs.GridSpec(ncols=10, nrows=13, figure=fig)
tscol = 7
seis     = fig.add_subplot(spec[:5,  :tscol])
mmseax   = fig.add_subplot(spec[5:10,:tscol])
crosssec = fig.add_subplot(spec[:10,    tscol:])
cb       = fig.add_subplot(spec[-1,  :tscol])

seis.plot(time, norm(trace), 'k')
# ================== END FORMAT OUTPUT FIGURE ======================


# For each slice we compute the MSE:
mMSE = []
print("Completed: ")
for i in range(no_slices):
    mmse_temp = ent.multiscale_entropy(time_series   = slices[i,:],
                                       sample_length = 2,
                                       tolerance     = tolerance,
                                       maxscale      = max_tau)
    mMSE.append(mmse_temp)
    # Plot MSE for slice
    crosssec.plot(m, mmse_temp, 'k', alpha=0.25)
    print(f"        {i+1}/{no_slices}")

# Format mMSE data
mMSE = np.array(mMSE)

# Output stats:
print("Max MMSE val:", np.max(mMSE.flatten()))
print("max value of any slice: ", np.nanmax(mMSE[np.logical_and(mMSE != np.nan, mMSE != np.inf) ]))

# Plot MMSE for all slices:
im = mmseax.imshow(X=np.transpose(mMSE),
             cmap='magma_r',
             norm=None,
             aspect='auto',
             extent=(time[0], time[-1], 1, max_tau),
             origin='lower',
             interpolation='none',
             vmin=0, vmax=vmax)

# Add colour bar
fig.colorbar(im, cax=cb, use_gridspec=True, orientation='horizontal')

# Fig bells and whistles:
#for axtmp in [mmseax, seis]:
#    axtmp.set_xlim([100, 1700])
seis.set_ylim([-1, 1])
mmseax.set_xlabel("Time [s]")
mmseax.set_ylabel("Tau")
seis.set_ylabel("Velocity")
crosssec.set_ylabel("SampEn")
crosssec.set_xlabel("Tau")
crosssec.set_ylim([0, 3])
cb.set_title(f'SampEn for {event}: dt = {dt}: tolerance = {tolerance} channel {channel} slice_time = {slice_time} s')
crosssec.set_ylim([0, 3])
seis.set_xlim(0,2000)
mmseax.set_xlim(0,2000)
print('setting x limits manually')

plt.show()
#plt.savefig(fname=f"./REV_MMSE_{event}_{channel}_{dt}_r{tolerance}_vmax{vmax}.pdf", dpi='figure', format='pdf')

