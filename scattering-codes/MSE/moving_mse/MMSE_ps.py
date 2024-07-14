from  obspy.core.trace import Trace
import matplotlib.gridspec as gs
from gen_MMSE_slices import gen_mmse_slices
import numpy as np
import pyentrp.entropy as ent
import matplotlib.pyplot as plt
from wetools import norm

# Input parameters:
event       = 'p-0.2_2hz_1_mfp_0.5_rad'   # event name
stn         = 't'                         # station
dt          =  '0.151'                    # downsampled timestep for time series
slice_time  = 50                          # seconds
max_tau     = 75                          # max scale
tolerance   = 0.03                        # mse tolerance
vmax        = 3.5                         # max value plot colourbar

# Data path
data_path   = f'../../../data/processed/p-0.2/{event}'

# Scale array
m = np.arange(1, max_tau+1)

# Number of timesteps in single slice
slice_length = int(slice_time/float(dt))

# Load and process data
event_np = np.loadtxt(f"{data_path}/M_0.ST{stn}.Z")
time_np  = np.loadtxt(f"{data_path}/time_data")

# Convert to np arrays:
tr = Trace()
tr.data = event_np
tr.stats.delta = time_np[1] - time_np[0]

# Convert to velocity
tr.differentiate()

# Resample at correct timestep
tr.resample(sampling_rate=(1 / float(dt)))
tr = tr.filter(type='bandpass', freqmin=0.1, freqmax=2)
trace = norm(tr.data)
length = len(trace)
time = np.linspace(0, float(dt)*length, length)


# Now need to generate slices for MMSE:
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
print("Completed: ")
for i in range(no_slices):
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

# Output stats:
print("max value of any slice: ", np.nanmax(mMSE[np.logical_and(mMSE != np.nan, mMSE != np.inf) ]))

# Plot MMSE for all slices:
im = mmseax.imshow(X=np.transpose(mMSE),
             cmap='magma_r',
             norm=None,
             aspect='auto',
             extent=(time[0], time[-1], 1, max_tau),
             origin='lower',
             vmin=0, vmax=vmax)

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
mmseax.set_xlim([0,300])
seis.set_xlim([0,300])

plt.show()
#plt.savefig(fname=f"./Figures/MMSE__PS.pdf", dpi='figure', format='pdf')

