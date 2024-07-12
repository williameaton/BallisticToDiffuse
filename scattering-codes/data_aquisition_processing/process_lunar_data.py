import matplotlib.pyplot as plt
import obspy
from wetools import obspy_gen_mpl, norm
# The .xml file comes from the Supp. Material of Nunn et al 2020:
# https://zenodo.org/record/3560482#.YriIGOzMIq1 - the DATALESS.zip file

# Set path to lunar raw data
path   = '../../data/lunar/'

# Define directories for the colloqial names used in paper
fullname = {'1971a': f"19710380047", '1971b': f"19712102100", '1972': f"19723452032"}
# Define channels
channelnames = {'lpx': 'MH1',
                'lpy': 'MH2',
                'lpz': 'MHZ',
                'spz': 'SHZ',
                }

# Define year (name) and station:
year = '1972'
stn  = '12'

# Flags:
filter_data = True         # Bandpass filter flag
plot_trace  = True         # Plot trace using MPL
cutoff_for_plot = 100000   # Timecutoff value for plotting - use large number to keep whole trace


# Set directory path:
folder = f"{path}/{fullname[year]}.gfs"
type = 'lpz'


# Read the lunar data:
data = obspy.read(f"{folder}/out.{fullname[year]}.{stn}.{type}.sac")
# Set channel name:
data[0].stats.channel = channelnames[type]
data[0].stats.network = 'XA'
data[0].stats.station = f'S{stn}'

# Convert to time series
x,y = obspy_gen_mpl(data[0])

# Remove insturment resp.
inv = obspy.read_inventory(f"{path}/responses/XA.1969-1977_updated_2019.xml",  format="STATIONXML")
data = data.remove_response(inventory=inv, output='VEL', pre_filt= (0.001, 0.005, 3, 4),
                    zero_mean=True, taper=True, taper_fraction=0.05, plot=True, fig=None)

#data.rotate(method='->ZNE', inventory=inv)
if filter_data:
    data = data.filter(type='bandpass', freqmin=0.1, freqmax=2)

X,Y = obspy_gen_mpl(data[0])

# Plot traces
if plot_trace:
    fig, ax = plt.subplots()
    ax.plot(x[x<cutoff_for_plot], norm(y[x<cutoff_for_plot]),'k')
    ax.plot(X[X<cutoff_for_plot], norm(Y[X<cutoff_for_plot]),'r', alpha=0.5)
    ax.legend('Before processing', 'After processing')
    plt.show()

if filter_data:
    filtstr = 'filtered_removed_response'
else:
    filtstr = 'NO_FILTER_removed_response'

#data.write(output_str = f"../../data/lunar/{filtstr}/{year}_{stn}_LP_{type}.mseed")
