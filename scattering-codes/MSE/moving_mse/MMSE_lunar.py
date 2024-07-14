import  obspy
import matplotlib.gridspec as gs
import numpy as np
import pyentrp.entropy as ent
import matplotlib.pyplot as plt
from gen_MMSE_slices import gen_mmse_slices
from wetools import norm

# Uncomment the event you want
event   = '1972'
stn     = '12'
t_start = 575
t_end   = 5375

"""
event   = '1971b'
stn     = '12'
t_start = 1780 - 100
t_end   = 6580 + 1000

event   = '1971a'
stn     = '14'
t_start = 485
t_end   = 5285
"""

# Set path to time series
data_path  = '../../../data/lunar/filtered_removed_response'

# Input parameters:
dt         = '0.151'   # downsampled timestep for time series
slice_time = 50        # window length in seconds
max_tau    = 100       # max scale
tolerance  = 0.03      # mse tolerance
vmax       = 3.5       # max value plot colourbar


# Create array of scales
m = np.arange(1, max_tau+1)

# Number of timesteps in single slice
slice_length = int(slice_time/float(dt))


# Load and process data
event_tr = obspy.read(f"{data_path}/{event}_{stn}_LP_lpz.mseed")[0]
event_tr.resample(sampling_rate=(1 / float(dt)))
event_tr = event_tr.slice(starttime=event_tr.stats.starttime, endtime=event_tr.stats.starttime+7000   )


# Convert to numpy arrays and normalise amplitude
length = len(event_tr.data)
time = np.linspace(0, float(dt)*length, length)
tmask = np.logical_and(time >= t_start - 200, time <= t_end)
time = time[tmask]
trace = norm(event_tr.data[tmask])

# Now need to generate slices for MMSE and time arrays:
timeslices, no_slices = gen_mmse_slices(trace   = time,
                                   slice_length = slice_length,
                                   overlap      = 0,
                                   normalise    = False)
slices, no_slices = gen_mmse_slices(trace       = trace,
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
print("completed: ")
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
print("max value of any slice: ", np.max(mMSE))

# Plot MMSE for all slices:
im = mmseax.imshow(X      =np.transpose(mMSE),
                   cmap   ='magma_r',
                   norm   = None,
                   aspect = 'auto',
                   extent = (time[0], time[-1], 1, max_tau),
                   origin = 'lower',
                   vmin   = 0,
                   vmax   = vmax)

# Add colour bar
fig.colorbar(im, cax=cb, use_gridspec=True, orientation='horizontal')


# Fig bells and whistles:
for axtmp in [mmseax, seis]:
    axtmp.set_xlim([time[0], time[-1]])
seis.set_ylim([-1, 1])
mmseax.set_xlabel("Time [s]")
mmseax.set_ylabel("Tau")
seis.set_ylabel("Velocity")
crosssec.set_ylabel("SampEn")
crosssec.set_xlabel("Tau")
cb.set_title(f'SampEn for {event}: dt = {dt}: tolerance = {tolerance} station {stn} slice_time = {slice_time} s')
crosssec.set_ylim([0, 3])

plt.show()
#plt.savefig(fname=f"./Figures/r{tolerance}_MMSE_{event}_{stn}_Z_{dt}_r{tolerance}_vmax{vmax}.pdf", dpi='figure', format='pdf')

