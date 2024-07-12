# Script to download Earth seismogram data from IRIS
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.signal.rotate import rotate2zne
from obspy.signal.filter import bandpass
import numpy as np

# Set Earth data parameters for download
output_name = 'earth_2013'
stime       = '2013-06-15T00:00:00.000'
etime       = '2013-06-15T01:00:00.000'

fmin = 0.1                          # Minimum frequency
fmax = 2                            # Maximum frequency

pre_filt = (0.005, 0.01, 2.0, 2.5)  # Prefilter before instrument response removal

print("Loading Earth data: ")
t1 = UTCDateTime(stime)
t2 = UTCDateTime(etime)

# Fetch waveform from IRIS FDSN web service into a ObsPy stream object
# and automatically attach correct response
fdsn_client = Client("IRIS")
st = fdsn_client.get_waveforms(network="HL", station="MHLO", location='--',
                               channel="HH?", starttime=t1, endtime=t2,
                               attach_response=True)


# define a filter band to prevent amplifying noise during the deconvolution
st.remove_response(output='VEL', pre_filt=pre_filt, plot=True)

# Rotate to NZE:
output_tuple = rotate2zne(data_1=st[0].data, azimuth_1=135.1, dip_1=-29.4,
                          data_2=st[1].data, azimuth_2=15, dip_2=-29.2,
                          data_3=st[2].data, azimuth_3=255, dip_3=-29.7)

# Assign BP-filtered, rotated data back to traces:
for i_update_index in range(3):
    # Bandpass the rotated data and re-assign to the stream
    st[i_update_index].data = bandpass(data=output_tuple[i_update_index],
                                       freqmin=fmin,
                                       freqmax=fmax,
                                       df=st[i_update_index].stats.sampling_rate,
                                       zerophase=True)

    # Uncomment if you want to upsample:
    # Up sample at 1000 Hz in accordance with parameter space data
    #st[i_update_index].interpolate(sampling_rate=samprate, method="linear")


print('Downloaded...\n', st)


# Writing to ascii and MSEED
np.savetxt(output_name + '.ascii', [st[0].data, st[1].data, st[2].data]) # 'int_S0173a_0.125_2BP.ascii'
print('Saved to ', output_name)
st.write(output_name + '.MSEED')
print('Written as MSEED: ', output_name + '.MSEED')

