from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.signal.rotate import rotate2zne
from obspy.signal.filter import bandpass

def download_mars(output_name, stime, etime, channel='BH?', fmin=0.125, fmax=2,
                  network='XB', station='ELYSE', location='02', client='IRIS',):

    print("Loading mars data: ")
    # Set time span for download
    t1 = UTCDateTime(stime)
    t2 = UTCDateTime(etime)

    # Fetch waveform from IRIS FDSN web service into a ObsPy stream object
    # and automatically attach correct response
    fdsn_client = Client(client)
    st = fdsn_client.get_waveforms(network=network, station=station, location=location,
                                   channel=channel, starttime=t1, endtime=t2,
                                   attach_response=True)

    print(f"sampling timestep: {st[0].stats.delta}")
    # define a filter band to prevent amplifying noise during the deconvolution
    st.remove_response(output='VEL',  pre_filt= (0.001, 0.005, 3, 4),
                       zero_mean=True, taper=True, taper_fraction=0.05)

    # Rotate using channel metadata accessed from: https://www.seis-insight.eu/en/science/seis-data/seis-metadata-access
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


    print('Downloaded...\n')
    print(st)

    st.write(output_filepath+output_name + f'_{fmin}_{fmax}.MSEED')
    print('Written as MSEED: ', output_filepath+output_name + f'_{fmin}_{fmax}.MSEED')



# Definitions for two different marsquake events
#output_name='S0173'
#stime='2019-05-23T02:15:00.000'
#etime='2019-05-23T03:15:00.000'

output_name = 'S0235'
stime = '2019-07-26T12:13:00.000'
etime = '2019-07-26T13:20:00.000'


# Uncomment for command line use
"""if __name__ == "__main__":
    # Get args from input:
    output_name = sys.argv[1]
    stime = sys.argv[2]
    etime = sys.argv[3]

    download_mars(output_name, stime, etime, channel='BH?', fmin=0.125, fmax=2, samprate=1000, network='XB',
                      station='ELYSE', location='02', client='IRIS', )
"""



output_filepath = '../../data/mars/filtered/'

download_mars(output_name, stime, etime, channel='BH?', fmin=0.1, fmax=2.0, network='XB',
                      station='ELYSE', location='02', client='IRIS' )