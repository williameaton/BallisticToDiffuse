import obspy
from calc_CC import calculate_CC
import numpy as np
import scipy.signal as ss
import matplotlib.pyplot as plt
from wetools import norm, obspy_gen_mpl

# Window length of 50 seconds
T = 50

# 500 s lapse time after start of trace
slicetimes = {'ARG':  {'2013': [22313, 22813], '2015': [432, 932] },
              'RDO':  {'2013': [703, 1203],    '2015': [494, 994] },
              'MHLO': {'2013': [700, 1200],    '2015': [205, 706] },
              }

# Load processed, rotated data:
data_path = f"../../data/earth/processed_rotated/"
station = 'ARG'

fmin = 0.1
fmax = 2

for event in ['2013', '2015']:
    print(f"{event} ")

    # Read in data and get number of channels:
    data = obspy.read(f"{data_path}/r{event}_{station}_vel.mseed")
    channels = len(data)

    fig, ax = plt.subplots(3,channels)

    # Loop each channel:
    for chl in range(channels):

        # Get trace data for single channel & convert to numpy array
        tr = data[chl]
        tr.data = norm(tr.data)
        print(f"  {tr.stats.channel} ")
        time, y = obspy_gen_mpl(tr)

        # Plot trace:
        ax[0,chl].plot(time, tr.data, 'k')

        # Slice the trace based on above-defined values and convert to numpy
        tr = tr.slice(starttime=tr.stats.starttime + slicetimes[station][event][0]  ,
                 endtime=tr.stats.starttime + slicetimes[station][event][1] )
        time, y = obspy_gen_mpl(tr)

        # Plot to emphasise sliced region
        ax[0, chl].plot(time + slicetimes[station][event][0], tr.data, 'r')

        # Compute CC of the sliced region
        freq_out, C, C_mean, no_windows = calculate_CC(tr.data, dt=tr.stats.delta, T=T)

        # ___________Calculate the power spectra for sliced time series___________

        # Taper data before FFT:
        taper = ss.tukey(M=len(tr.data), alpha=0.1, sym=False)
        tracetaper = taper * tr.data
        fft = np.abs(np.fft.fft(tracetaper, n=len(tracetaper) * 2))
        fftfreq = np.fft.fftfreq(n=len(tracetaper) * 2, d=float(tr.stats.delta))

        # Define mask of frequencies between fmin and fmax Hz
        m = np.logical_and(fftfreq>=fmin, fftfreq<=fmax)

        # Uncomment to save FFT data:
        # np.savetxt(fname=f'NEW_mars_power_LP_S0173_chl{chls[i]}_delta_{np.around(downsample_dt,3)}',
        #           X=[fftfreq, fft])

        # Plot spectrum
        ax[2,chl].plot(fftfreq[m], fft[m])
        ax[2,chl].set_xlim([fmin, fmax])
        #_______________________________________________________


        # Slice out frequencies between fmin and fmax
        indl = np.where(freq_out>=fmin)[0][0]
        indh = np.where(freq_out > fmax)[0][0]
        freq_out = freq_out[indl:indh]
        C = C[indl:indh, indl:indh]
        cm = np.mean(C[np.triu_indices_from(C, k=1)])
        print("cm: ", cm)

        # Plot matrix:
        im = ax[1,chl].imshow(C,
                              extent=(freq_out[0], freq_out[-1], freq_out[0], freq_out[-1]),
                              origin='lower',
                              cmap='magma_r')

        # Set plot title and colourbar
        fig.colorbar(mappable=im, ax=ax[1,chl], use_gridspec=True)
        ax[0,chl].set_title(f"channel  {tr.stats.channel}  Mean CC: {cm}", fontsize=4)

    # Save figure
    plt.savefig(fname=f'./figs/cc_{event}_{station}_T{T}.pdf', dpi='figure', format='pdf')