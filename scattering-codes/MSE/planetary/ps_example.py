# Example MSE calculation for PS data
from  obspy.core.trace import Trace
import scipy.signal as ss
import numpy as np
import pyentrp.entropy as ent
import scipy.integrate as spi
import matplotlib.pyplot as plt
from wetools import norm


# Input parameters:
event     = 'p-0.2_2hz_1_mfp_0.5_rad'
stn       = 't'                                                         # select single station
data_path = f'../../../../Scattering/data/processed/p-0.2/{event}'
chl       = 'Z'
dt        = '0.151'                                                     # downsampled timestep
max_tau   = 200                                                         # Max scale
r         = 0.01                                                        # threshold
m = np.arange(1, max_tau+1)                                             # Scales


# Load and process data
event_np = np.loadtxt(f"{data_path}/M_0.ST{stn}.{chl}")
time_np  = np.loadtxt(f"{data_path}/time_data")

# Load data to trace
tr = Trace()
tr.data = event_np
tr.stats.delta = time_np[1] - time_np[0]

# Convert to velocity
tr.differentiate()

# Resample
tr.resample(sampling_rate=(1 / float(dt)))
trace = norm(tr.data)

# Create time array for downsampled array
length = len(trace)
time = np.linspace(0, float(dt)*length, length)

# Compute power spectrum for reference
taper = ss.tukey(M=len(trace), alpha=0.1, sym=False)
tracetaper = taper*trace
fft = np.abs(np.fft.fft(tracetaper, n=len(tracetaper)*2))
fftfreq = np.fft.fftfreq(n=len(tracetaper)*2, d=float(dt))
np.savetxt(fname=f'ps_power',
           X=[fftfreq, fft])

# Compute MSE and save
MSE = ent.multiscale_entropy(time_series=trace, sample_length=2, tolerance=r, maxscale=200)
np.savetxt(fname=f'PS _{event}_r_{r}_example_{np.around(float(dt),3)}' , X=MSE)

# Calculate integrated MSE
Tm = spi.trapz(MSE, m)
print(f"Tm:  {Tm}")


figMSE, axMSE = plt.subplots()
axMSE.plot(m, MSE)
axMSE.set_title(f"Integrated m: {Tm}")
axMSE.legend([f'Channel: {chl}'])

plt.show()

