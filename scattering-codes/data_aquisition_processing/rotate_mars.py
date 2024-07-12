# Rotates processed Mars data stored in `data/mars/filtered`
# to ZNE and stores in RTZ in data/processed_rotated
import matplotlib.pyplot as plt
import obspy
from rotate import rotate_stream
import numpy as np

# Set data path and frequency limits of trace you want to load
data_path = f"../../data/mars/filtered"
fmin = 0.1
fmax = 0.6
# Set InSight station coordinates for rotation
stn_coords = [4.5024, 135.6234]

# Uncomment for the event you want to use:
event = 'S0173_NEW'
ev_lat = 3.398
ev_lon = 165
ev_az  = 91
ev_dis = 30

"""event = 'S0235_NEW'
ev_lat = 11.5
ev_lon = 163
ev_az  = 74
ev_dis = 28.75"""


data = obspy.read(f'{data_path}/{event}_{fmin}_{fmax}.MSEED')

# Currently data are downloaded in ZNE but labelled UVW - assuming that is equivalent to ZNE now rotated:
data.select(channel='BHV')[0].stats.channel = 'N'
data.select(channel='BHW')[0].stats.channel = 'E'

s = data[0].stats
x = np.linspace(0, s.npts*s.delta, s.npts)


# Plot before and after rotation
fig, ax = plt.subplots(3)
for axi in range(3):
    ax[axi].plot(x[x>200], data[axi].data[x>200], 'k')

rotate_stream(stream=data,
              method="NE->TP",
              src=[ev_lat, ev_lon],
              stn=stn_coords,
              geoco=1-0.00589,
              overwrite_channel_names=False, invert_p=False, invert_t=False)

# Update channel names?
data.select(channel='BHU')[0].stats.channel = 'Z'
data.select(channel='N')[0].stats.channel   = 'R'
data.select(channel='E')[0].stats.channel   = 'T'

for axi in range(3):
    ax[axi].plot(x[x>200], data[axi].data[x>200], 'r')
    ax[axi].set_title(f"Channel: {data[axi].stats.channel}")
plt.show()

# Write rotated
data.write(f"{data_path}/{event}r_{fmin}_{fmax}.MSEED")
