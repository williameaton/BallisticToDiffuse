import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append("../")
import obspy as ob
import obspy.core.stream as ocs
from load_avg_Qc import calc_avg_QC
from qc_funcs import plot_decay_curve, get_QC

# SPECIFY EVENT and time slices used to compute Qc:
# Uncomment the desired event
# Need to use these slices for this one because of glitches
"""
event = "1971b"
stn   = "12"
slice_start = 1780
slice_end   = 6580
offset = 220
"""
"""
event = "1971a"
stn   = "14"
slice_start = 485
slice_end   = 5285
offset = 500

"""
event = "1972"
stn   = "12"
slice_start = 575
slice_end   = 5375
offset = 210

period = 'LP'

# Set data path
data_path = f"../../../data/lunar/filtered_removed_response/"

# Read 3 channels
stream  = ocs.Stream()
for chl in ["x","y","z"]:
    fname = f"{event}_{stn}_{period}_lp{chl}.mseed"
    stream += ob.read(f"{data_path}/{fname}")[0]

# Get stream stats
stats = stream[0].stats

# Slice the streams between slice_start and slice_end
stream.slice(starttime=stats.starttime + slice_start, endtime=stats.starttime + slice_end)
time = np.linspace(0, stats.npts*stats.delta, stats.npts)

# Compute Qc and plot figure
fig, ax, Qc_arr, time, energy = get_QC(time=time,
                                        R=stream[0].data,
                                        T=stream[1].data,
                                        Z=stream[2].data,
                                        convert_to_vel=False)

# Array index for max energy (normalised)
ind = np.where(energy == 1)[0][0]
mean_orig = np.mean(Qc_arr[:, 2])
std_orig = np.std(Qc_arr[:, 2])

# Plot original decay:
ax[2] = plot_decay_curve(ax=ax[2], mean=mean_orig, std=std_orig, time=time, ind=ind, freq=1,
                         offset=offset)


# Calculate reduced mean/std
mean_reduced, std_reduced = calc_avg_QC(d=Qc_arr, upper_cutoff=10000)
ax[3] = plot_decay_curve(ax=ax[3], mean=mean_reduced, std=std_reduced, time=time, ind=ind,
                         freq=1, offset=offset)

# Set titles
ax[0].set_title(f"{event}  {stn}")
ax[2].set_title(f"{mean_orig} +- {std_orig}: offset: {offset}")
ax[3].set_title(f"{mean_reduced} +- {std_reduced}: offset: {offset}")


fig.suptitle(f"Event: {event}   Station: {stn}")
for i in range(4):
    ax[i].set_ylabel("Norm. Energy")

ax[-1].set_xlabel("Time [s]")

print(mean_reduced, std_reduced)


plt.savefig(fname=f"./pdfs/REV_LUNAR_QC_{event}_{stn}.pdf", dpi='figure', format='pdf')
#plt.show()