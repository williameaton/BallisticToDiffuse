# Writing these individually for the Mars data because I want to avoid glitches
import obspy
import obspy.core.trace as oct
import obspy.core.stream as ocs
import matplotlib.pyplot as plt
import scipy.signal as ss
import numpy as np
import scipy.integrate as spi
from wetools import norm
import pyentrp.entropy as ent

# For the lunar and Martian MSE we want to look at the effect of
# band-passing the trace before computing the MSE.
# here we loop over a number of different frequency bands and
# compute the MSE for each

# Define frequency minimum and maximums:
F = [[0.1, 2.0]]#, [0.1, 0.6], [0.6, 1.1], [1.1, 1.75], [1.75, 2.0]]

# MSE stuff:
r = 0.03                        # MSE threshold
m = np.arange(1, 201)           # Scale
downsample_dt = 0.151           # Resample timestep (downsample data) -- long period lunar
chls = ['Z', 'R', 'T']


starttime = 100
# Defining start/end of three sections separated by glitches:
# Offset used to pick it manually
s1 = 270
e1 = 1342
s2 = 1378
e2 = 1645.6
s3 = 1661
e3 = 2271.4
times = [[s1, e1], [s2, e2], [s3, e3]]

for MM in range(len(F)):
    # Define bandpass
    fmin = F[MM][0]
    fmax = F[MM][1]
    print(f"------------ FREQ BAND: {fmin} to {fmax} --------------")

    # Load data - has already been filtered
    file_path = f'../../../data/mars/filtered/S0235_NEWr_{fmin}_{fmax}.MSEED'
    st = obspy.read(file_path)

    # Original slicing of data
    st = st.slice(starttime=st[0].stats.starttime+starttime, endtime=st[0].stats.starttime+2502)

    # Initialise arrays to store data
    glitched   = []
    deglitched = []

    for chl in chls:
        print(f'Channel: {chl}')

        # Used to concatenate deglitched data
        new_trace = np.empty(shape=(1))*0

        # Get data for specific channe;
        tr = st.select(channel=chl)[0]

        # Normalise trace by maximum amplitude and get dt and
        # number of samples in time series
        trace = norm(tr.data)
        stats = tr.stats
        dt = stats.delta
        n  = stats.npts

        # Generate time array corresponding to trace
        time = np.linspace(0, n*dt, n)

        # Mask of data between beginning of section 1 and end of section 3
        # Add it to 'glitched' array
        mask = np.logical_and(time>s1, time<e2)
        glitched.append(norm(trace[mask]))

        # Look though sections to concatenate parts of the data without glitches
        # and store in new_trace
        # Note that this does NOT leave time gaps where the glitches were
        for j in range(len(times)):
            mask = np.logical_and(time>times[j][0], time<times[j][1])
            new_trace = np.concatenate((new_trace, trace[mask]))

        # create deglitched new time array and plot :
        new_time = np.linspace(0, len(new_trace)*dt, len(new_trace)) + s1

        # Store deglitched
        # Needs to be re-normed post deglitching
        new_trace = norm(new_trace)
        deglitched.append(new_trace)
        # ---------------- END OF LOOP OVER CHANNELS ----------------

    # Now we have all the time series containing glitches (glitched)
    # and with them cut out (deglitched)
    # convert to numpy
    glitched = np.array(glitched)
    deglitched = np.array(deglitched)


    # Downsample the de-glitched data for MSE processing
    # --> first, add it to obspy stream
    deglitched_st = ocs.Stream()
    for i in range(3):
        deglitched_tr = oct.Trace()
        deglitched_tr.data = deglitched[i,:]
        deglitched_tr.stats.delta = dt
        deglitched_tr.stats.channel = chls[i]
        deglitched_st += deglitched_tr
    # --> now run the downsampling
    deglitched_st.resample(sampling_rate=(1/downsample_dt))

    # Create new array of down-sampled data
    downsampled = []
    for i in range(len(chls)):
        downsampled.append(deglitched_st[i].data)
    downsampled = np.array(downsampled)

    # Save downsampled, deglitched data for for MMSE analysis
    np.savetxt(fname=f'../../../data/mars/deglitched_downsampled/HREV_NEW_S0235_delta_{np.around(downsample_dt,3)}_{fmin}_{fmax}', X=downsampled)

    # Create time array for the downsampled, deglitched data
    # starts at time s1
    downsample_time = np.linspace(0, len(downsampled[0,:])*downsample_dt, len(downsampled[0,:]) ) + s1


    # Now calculate MSE:
    print("Running MSE ")
    for i in range(3):

        # Find the time of max energy within the downsampled,deglitched time series
        tmax = downsample_time[np.where(downsampled[i,:] == np.max(downsampled[i,:]))[0][0]]

        # Take some time frame of data after the max amplitude arrival
        tmask = np.logical_and(downsample_time >= tmax, downsample_time<= tmax +1500)

        # Select this part of downsampled timeseries
        trace = downsampled[i,tmask]

        # Compute frequency content of trace for reference (later plotted in Fig 10. b)
        for j in range(len(chls)):
            taper = ss.tukey(M=len(trace), alpha=0.1, sym=False)
            tracetaper = taper * trace
            fft = np.abs(np.fft.fft(tracetaper, n=len(tracetaper) * 2))
            fftfreq = np.fft.fftfreq(n=len(tracetaper) * 2, d=float(dt))
            np.savetxt(fname=f'./mars/REV_mars_power_LP_S0235_chl{chls[j]}_delta_{np.around(downsample_dt,3)}',
                       X=[fftfreq, fft])

        # Compute MSE
        MSE = ent.multiscale_entropy(time_series=norm(trace), sample_length=2, tolerance=r, maxscale=200)

        # Save MSE TRACE FOR FREQ COMPARISON:
        np.savetxt(f'./mars/JREV_FINAL_MSE_S0235_{chls[i]}_{fmin}_{fmax}', X=MSE)

        # Calculate the area under m-curve (intMSE):
        Tm = spi.trapz(MSE, m)
        print(f"Tm:  {Tm}")

        # Redundant save
        #np.savetxt(fname=f'./mars/WE_LP_S0235_chl{chls[i]}_delta_{np.around(downsample_dt,3)}_r{r}' , X=MSE)
