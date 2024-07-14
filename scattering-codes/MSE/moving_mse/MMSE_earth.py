import os
import obspy
import matplotlib.gridspec as gs
import numpy as np
import pyentrp.entropy as ent
import matplotlib.pyplot as plt
from gen_MMSE_slices import gen_mmse_slices
from wetools import norm

# Start times of slices
slicetimes = {'ARG': {'2013' : [22313], '2015': [432]},
              'RDO': {'2013' : [703],   '2015': [494]},
              'MHLO': {'2013': [700],   '2015': [206]}
              }

# Set path to time series
data_path = '../../../data/earth/processed_rotated/'

# Input parameters:
event      = '2015'     # event year
station    = 'MHLO'     # station name
dt         = '0.151'    # downsampled timestep for time series
channel    = 'Z'        # data channel

slice_time = 50         # window length in seconds
max_tau    = 100        # max scale
tolerance  = 0.03       # mse tolerance
vmax       = 3.5        # max value plot colourbar


# Create array of scales
m = np.arange(1, max_tau+1)

# Define start and tend time of trace
stime = slicetimes[station][event][0] - 100
etime = stime + 800

# Number of timesteps in single slice
slice_length = int(slice_time/float(dt))

# Load and process data
event_tr = obspy.read(f"{data_path}/r{event}_{station}_vel.MSEED").select(channel=channel)[0]
event_tr.resample(sampling_rate=(1 / float(dt)))
event_tr = event_tr.slice(starttime=event_tr.stats.starttime +stime ,
                          endtime=event_tr.stats.starttime + etime)

# Convert to numpy arrays and normalise amplitude
trace = norm(event_tr.data)
length = len(trace)
time = np.linspace(0, float(dt)*length, length)


# Now need to generate slices for MMSE and time arrays:
timeslices, no_slices = gen_mmse_slices(trace   = time,
                                   slice_length = slice_length,
                                   overlap      = 0,
                                   normalise    = False)
slices, no_slices      = gen_mmse_slices(trace  = trace,
                                   slice_length = slice_length,
                                   overlap      = 0)


# ==================== FORMAT OUTPUT FIGURE: ======================
fig = plt.figure(constrained_layout=True, figsize=(13,6.5))
spec = gs.GridSpec(ncols=10, nrows=13, figure=fig)
tscol = 7
seis     = fig.add_subplot(spec[:5,  :tscol])
mmseax   = fig.add_subplot(spec[5:10,:tscol])
crosssec = fig.add_subplot(spec[:10,    tscol:])
cb       = fig.add_subplot(spec[-1,  :tscol])

seis.plot(time, trace, 'k')
# ================== END FORMAT OUTPUT FIGURE ======================


# For each slice we compute the MSE:
mMSE = []

# Loop each slice
print("completed ...")
for i in range(no_slices):
    # Compute mse for slice and store in mMSE array
    mmse_temp = ent.multiscale_entropy(time_series=slices[i,:],
                                       sample_length=2,
                                       tolerance=tolerance,
                                       maxscale=max_tau)
    mMSE.append(mmse_temp)
    # Plot MSE for slice
    crosssec.plot(m, mmse_temp, 'k', alpha=0.25)
    print(f"        {i+1}/{no_slices}")

# Format mMSE data
mMSE = np.array(mMSE)
mmm = mMSE.flatten()

# Output stats:
print("Max MMSE val:", np.max(mmm))
print("Max MMSE val without nan or inf:", np.max(mmm[np.logical_and(mmm<np.inf,~np.isnan(mmm))]))
print("max value of any slice: ", np.nanmax(mMSE[np.logical_and(mMSE != np.nan, mMSE != np.inf) ]))


# Plot MMSE for all slices:
im = mmseax.imshow(X      = np.transpose(mMSE),
                   cmap   = 'magma_r',
                   norm   = None,
                   aspect = 'auto',
                   extent = (time[0], time[-1], 1, max_tau),
                   origin = 'lower',
                   vmin   = 0,
                   vmax   = vmax)

# Add colour bar
fig.colorbar(im, cax=cb, use_gridspec=True, orientation='horizontal')

# Figure bells and whistles:
for axtmp in [mmseax, seis]:
    axtmp.set_xlim([time[0], time[-1]])
seis.set_ylim([-1, 1])
mmseax.set_xlabel("Time [s]")
mmseax.set_ylabel("Tau")
seis.set_ylabel("Velocity")
crosssec.set_ylabel("SampEn")
crosssec.set_ylim([0, 3])
crosssec.set_xlabel("Tau")
cb.set_title(f'SampEn for {event}: dt = {dt}: tolerance = {tolerance} channel {channel} slice_time = {slice_time} s')

crosssec.set_ylim([0, 3])

strn = f"./Figures/r{tolerance}/MMSE_{event}_{station}_{channel}_{dt}_r{tolerance}_vmax{vmax}.pdf"

# Save figure or make directory if necessary
try:
    plt.savefig(fname=strn, dpi='figure', format='pdf')
except:
    os.mkdir(f"./Figures/r{tolerance}")
    plt.savefig(fname=strn, dpi='figure', format='pdf')

plt.show()