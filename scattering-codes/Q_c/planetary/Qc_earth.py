import obspy
import numpy as np
import sys; sys.path.append('../')
from calc_qc_func import get_QC
from qc_funcs import plot_decay_curve
from load_avg_Qc import calc_avg_QC

networks = {"ARG": "HL", "RDO": "HL", "ITM": "HL", "VSU": "GE", "A36M": "TA", "MHLO": "HL"}

# Time frame for timeseries in which Qc is being calculated
slicetimes = {'ARG':  {'2013': [22313, 22813], '2015': [432, 932]},
              'RDO':  {'2013': [703,   1203],  '2015': [494, 994]},
              'MHLO': {'2013': [700,   1200],  '2015': [205, 706]},
              'SIVA': {'2013': [676,   1176],  '2015': [195, 695]}
              }

data_path = '../../../data/earth'    # specify data path
stn       = 'ARG'                    # station name
offset    = 5                        # offset in time just for plotting


# Loop different events specified by their year
for year in ['2015']:

    # Load processed Earth data
    st = obspy.read(f"{data_path}/processed_rotated/r{year}_{stn}_vel.mseed")

    # Trim data between specified times (slicetimes)
    stats = st[0].stats
    st = st.trim(starttime=stats.starttime + slicetimes[stn][year][0],
                 endtime=stats.starttime   + slicetimes[stn][year][1])

    # Sanity check
    assert(st[0].stats.delta == st[1].stats.delta )
    assert(st[0].stats.delta == st[2].stats.delta )
    assert(st[0].stats.npts  == st[1].stats.npts  )
    assert(st[0].stats.npts  == st[2].stats.npts  )

    # Create time array
    time = np.linspace(0, stats.npts*stats.delta, stats.npts)

    # Calculating energy:
    fig_result, ax_result, D, time, energy = get_QC(time, st[0].data, st[1].data, st[2].data, convert_to_vel=False)

    # Array index for max energy (normalised)
    ind = np.where(energy == 1)[0][0]
    mean_orig = np.mean(D[:, 2])
    std_orig = np.std(D[:, 2])


    # Plot original decay:
    ax_result[2] = plot_decay_curve(ax=ax_result[2], mean=mean_orig, std=std_orig, time=time, ind=ind, freq=1, offset=offset)


    # Calculate reduced mean/std
    mean_reduced, std_reduced = calc_avg_QC(d=D, upper_cutoff=400)
    ax_result[3] = plot_decay_curve(ax=ax_result[3], mean=mean_reduced, std=std_reduced, time=time, ind=ind, freq=1, offset=offset)

    # Set titles
    ax_result[0].set_title(f"{stn}  {year}")
    ax_result[2].set_title(f"{mean_orig} +- {std_orig}: offset: {offset}")
    ax_result[3].set_title(f"{mean_reduced} +- {std_reduced}: offset: {offset}")

    # uncomment to show or to save as a PDF
    #plt.show()
    #plt.savefig(fname=f"./pdfs/EARTH_QC_{stn}_{year}.pdf", dpi='figure', format='pdf')


