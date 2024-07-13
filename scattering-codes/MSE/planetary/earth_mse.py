import obspy
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi
import scipy.signal as ss
from wetools import norm, obspy_gen_mpl
import pyentrp.entropy as ent

# MSE stuff:
r = 0.03                    # Threshold value
m = np.arange(1,201)        # Scale values

# Earth ALWAYS NEEDS DOWNSAMPLING to be comparable with lunar
downsample_dt = 0.151 # long period lunar

# Define station and event
station = 'ARG'
event = '2013'

# Set data path and channels to loop through
chls = ['Z']
file_path = '../../../data/earth/processed_rotated/'

# --------------------------------------------------------------------------------------------------------

lapse = 500
# Defining slice start time
slicetimes = {'ARG': {'2013' : [22313], '2015': [432]},
              'RDO': {'2013' : [703],   '2015': [494]},
              'MHLO': {'2013': [700],   '2015': [206]},
              'SIVA': {'2013': [676, 1176], '2015': [195, 695]}
              }

# Load data and slice
st = obspy.read(f"{file_path}/r{event}_{station}_vel.mseed")
st_original = st.copy()
st = st.slice(starttime=st[0].stats.starttime + slicetimes[station][event][0],
              endtime=st[0].stats.starttime   + slicetimes[station][event][0]+ 10000
              )

# Normalise the data:
for i in range(3):
    st[i].data = norm(st[i].data)

# Loop channels
for chl in chls:
    # Get trace for specific channel & convert to numpy
    tr_original = st_original.select(component=chl)[0]
    tr = st.select(component=chl)[0]
    x, y = obspy_gen_mpl(tr)

# Now we have sliced data we need to downsample the trace
# by first putting it into a stream
st.resample(sampling_rate=(1/downsample_dt))
downsampled = []
for chl in chls:
    downsampled.append(norm(st.select(channel=chl)[0].data))
# Convert back to numpy
downsampled = np.array(downsampled)


# Now calculate MSE:
print("Running MSE ")

# tms stores Tm values for all channels
tms=[]
for i in range(len(chls)):

    # Take only stuff after the maximum amplitude:
    down_y = norm(downsampled[i,:])
    # Length of downsampled tdata
    lendown = len(down_y)
    # Create time array for downsampled data
    downsampled_time = np.linspace(0, lendown*downsample_dt, lendown) + slicetimes[station][event][0]
    # Get location of max amplitude
    maxamploc = downsampled_time[np.where(np.abs(down_y)==np.max(np.abs(down_y)))[0][0]]

    # Mask of times between max amplitude and max amplitude time + lapse time
    # Plot this
    timemask = np.logical_and(downsampled_time>= maxamploc, downsampled_time <= maxamploc+lapse)
    xplot = downsampled_time[timemask]
    yplot = down_y[timemask]

    # Compute frequency content of trace for reference (later plotted in Fig 10. b)
    taper = ss.tukey(M=len(yplot), alpha=0.1, sym=False)
    tracetaper = taper * yplot
    fft = np.abs(np.fft.fft(tracetaper, n=len(tracetaper) * 2))
    fftfreq = np.fft.fftfreq(n=len(tracetaper) * 2, d=float(downsample_dt))
    # Save frequency content to file:
    np.savetxt(fname=f'./earth/REV_earth_power_{station}{event}_chl{chls[i]}_0.1_2.0',
               X=[fftfreq, fft])

    # Compute the MSE
    MSE = ent.multiscale_entropy(time_series=yplot, sample_length=2, tolerance=r, maxscale=200)

    # Calculate the area under m-curve
    # (integrated MSE)
    Tm = spi.trapz(MSE, m)
    # Store integrated MSE
    tms.append(Tm)

    # Save MSE array
    outstr = f'./earth/REV_LP_{station}{event}_chl{chls[i]}_delta_{np.around(downsample_dt,3)}_r{r}'
    np.savetxt(fname=outstr, X=MSE)
    print(f'Saved MSE to {outstr}')



# Tm holds integrated MSE -- save here if you want?
print(f"Integrated MSE values: ", tms)

# Output MSE:
print(f"max mse: {np.nanmax(MSE)} at {m[np.where(MSE==np.nanmax(MSE))]}")

