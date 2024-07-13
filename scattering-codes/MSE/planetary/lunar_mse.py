import obspy
import numpy as np
import scipy.integrate as spi
import scipy.signal as ss
import pyentrp.entropy as ent
from wetools import norm

# For the lunar and Martian MSE we want to look at the effect of
# band-passing the trace before computing the MSE.
# here we loop over a number of different frequency bands and
# compute the MSE for each

# Define frequency minimum and maximums:
F = [[0.1, 2] , [0.1, 0.6], [0.6, 1.1], [1.1, 1.75], [1.75, 2.0]]

# Input parameters:
data_path = '../../../data/lunar/NO_FILTER_removed_response/'
dt = '0.151'                    # Resample timestep (downsample data)
r = 0.03                        # MSE threshold
chl = 'z'                       # Desired channel (file name)   TODO: combine this with below?
CHL = 'MHZ'                     # Obspy channel name
m = np.arange(1, 201)           # Scale

# Uncomment the event you want to use
"""event = '1972'
stn      = '12'
t_start  = 575
t_end    = 10000"""

"""event = '1971b'
stn      = '12'
t_start  = 1780
t_end    = 10000"""

event    = '1971a'
stn      = '14'
t_start  = 485
t_end    = 7200



# Loop through frequency bands
for MM in range(len(F)):
    # Define bandpass
    fmin = F[MM][0]
    fmax = F[MM][1]
    print(f"------------ FREQ BAND: {fmin} to {fmax} --------------")

    # Load data + filter and resample
    event_tr = obspy.read(f"{data_path}/{event}_{stn}_LP_lp{chl}.mseed").select(channel=CHL)[0]
    event_tr = event_tr.filter(type='bandpass', freqmin=fmin, freqmax=fmax)
    event_tr.resample(sampling_rate=(1 / float(dt)))

    # Convert to numpy and generate time arrays
    trace = norm(event_tr.data)
    length = len(trace)
    time = np.linspace(0, float(dt)*length, length)

    # slice trace:
    print(f"Slicing from {t_start} to {t_end}")
    mask = np.logical_and(time > t_start, time < t_end)
    trace = trace[mask]
    time  = time[mask]

    # Fetch only the stuff after the maximum amplitude:
    mintime = time[np.where(np.abs(trace)==np.max(np.abs(trace)))[0][0]]
    tmask = np.logical_and(time >= mintime, time<= mintime + 4500)
    print(f"max amp at {mintime}")

    # Compute frequency content of trace for reference (later plotted in Fig 10. b)
    taper = ss.tukey(M=len(trace[tmask]), alpha=0.1, sym=False)
    tracetaper = taper*trace[tmask]
    fft = np.abs(np.fft.fft(tracetaper, n=len(tracetaper)*2))
    fftfreq = np.fft.fftfreq(n=len(tracetaper)*2, d=float(dt))
    np.savetxt(fname=f'./lunar/REV_lunar_{chl}_power_{stn}{event}_chl{chl}_{fmin}_{fmax}',
               X=[fftfreq, fft])

    # Compute MSE
    MSE = ent.multiscale_entropy(time_series=norm(trace[tmask]), sample_length=2, tolerance=r, maxscale=200)

    # Calculate the area under m-curve (intMSE):
    Tm = spi.trapz(MSE, m)

    # Save MSE data to file
    outstr = f'./lunar/REV_LP_{stn}{event}_chl{chl}_delta_{np.around(float(dt),3)}_r{r}'
    np.savetxt(fname=outstr , X=MSE)
    print(f'Saved MSE to {outstr}')

    # Save MSE TRACE FOR FREQ COMPARISON:
    # Second save is redundant but stores the frequency bands
    np.savetxt(f'./lunar/REV_MSE_{event}_{stn}_Z_{fmin}_{fmax}', X=MSE)

