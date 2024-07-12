# Rotates processed Earth data stored in data/earth to ZNE and stores in RTZ in data/processed_rotated
import matplotlib.pyplot as plt
import obspy
from rotate import rotate_stream
import numpy as np

# Define station coordinates
stn_coords = {'ARG':  [36.213558, 28.12122],
              'RDO':  [41.145031, 25.33553],
              'MHLO': [36.68984,  24.40171],
              }

# Define event coordinates
event_coords = {'2013': [34.22, 25.06],
                '2015': [35.2337, 26.82],
                '2020': [37.918, 26.790]}

data_path = f"../../data/earth/"
year = '2013'
stn  = 'ARG'

# Read in mseed data
data = obspy.read(f'{data_path}/processed/{year}_{stn}_vel_newBP.mseed')

# Currently data are downloaded in ZNE
fig, ax = plt.subplots(3)
s = data[0].stats
x = np.linspace(0, s.npts*s.delta, s.npts)

for axi in range(3):
    ax[axi].plot(x, data[axi].data, 'k')


rotate_stream(stream=data,
              method="NE->TP",
              src=event_coords[year],
              stn=stn_coords[stn],
              geoco=1-0.00335,
              overwrite_channel_names=False, invert_p=False, invert_t=False)

data.select(channel='HHZ')[0].stats.channel = 'Z'
data.select(channel='HHN')[0].stats.channel = 'R'
data.select(channel='HHE')[0].stats.channel = 'T'

for axi in range(3):
    ax[axi].plot(x, data[axi].data, 'r')
    ax[axi].set_title(f"{data[axi].stats.channel}")
plt.show()

# Write rotated
data.write(f"{data_path}/processed_rotated/r{year}_{stn}_vel.mseed")
