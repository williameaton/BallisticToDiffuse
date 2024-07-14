import numpy as np
import matplotlib.pyplot as plt
import obspy
from wetools import norm
import sys; sys.path.append('../')
from load_avg_Qc import calc_avg_QC
from qc_funcs import plot_decay_curve, get_QC_interpolate

def calc_exp_curve(t_exp, freq, Qc):
    return norm(np.exp(-(2 * np.pi * freq * t_exp) / Qc))

# Load data that still has glitches. We then set those glitch values = 0 rather than
# cutting them out to preserve time decay

# Filtered data path
data_path = f"../../../data/mars/filtered"

event_name = 'S0173'

if event_name=='S0173':
    offset = 100
    starttime = 600

    # Define glitch times:
    s1 = 1366 - starttime
    e1 = 1417.5 - starttime
    s2 = 1814 - starttime
    e2 = 1839.5 - starttime
    glitches = [[s1, e1], [s2,e2]]
elif event_name=='S0235':
    offset = 0
    starttime = 100
    s1 = 682
    e1 = 695
    s2 = 1342
    e2 = 1378
    s3 = 1645.6
    e3 = 1661
    glitches = [[s1, e1], [s2,e2], [s3,e3]]
else:
    raise ValueError('Must be S0173 or S0235')


cutoff = 2000
dt     = 0.05

# Read in data and slice the correct time length
d = obspy.read(f'{data_path}/{event_name}_NEWr_0.1_2.0.MSEED')
d = d.slice(starttime=d[0].stats.starttime + starttime, endtime=d[0].stats.starttime+2000)

# Select channel data
Z = d.select(channel='Z')[0]
time = np.linspace(0, Z.stats.delta*Z.stats.npts, Z.stats.npts) + starttime

R = d.select(channel='R')[0].data
T = d.select(channel='T')[0].data
Z = Z.data


# Interpolate the values to 0 for glitches and compute Qc
fig_result, ax_result, D, time, energy  = get_QC_interpolate(time, R, T, Z,
                                                             convert_to_vel=False,
                                                             glitches=glitches,
                                                             starttime=starttime)

# Array index for max energy (normalised)
ind = np.where(energy == 1)[0][0]
mean_orig = np.mean(D[:, 2])
std_orig = np.std(D[:, 2])


# Plot original decay:
ax_result[2] = plot_decay_curve(ax=ax_result[2],
                                mean=mean_orig,
                                std=std_orig,
                                time=time,
                                ind=ind,
                                freq=1,
                                offset=offset)


# Calculate reduced mean/std
mean_reduced, std_reduced = calc_avg_QC(d=D, upper_cutoff=cutoff)
ax_result[3] = plot_decay_curve(ax=ax_result[3],
                                mean=mean_reduced,
                                std=std_reduced,
                                time=time,
                                ind=ind,
                                freq=1,
                                offset=offset)

# Set figure titles
ax_result[0].set_title(f"{event_name}")
ax_result[2].set_title(f"{mean_orig} +- {std_orig}: offset: {offset}")
ax_result[3].set_title(f"{mean_reduced} +- {std_reduced}: offset: {offset}")


#plt.show()
plt.savefig(fname=f"./pdfs/MARS_QC_{event_name}_GLITCHES_INTERP.pdf", dpi='figure', format='pdf')
