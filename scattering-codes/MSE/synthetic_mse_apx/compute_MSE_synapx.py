# Computes MSE for synthetic time series described in Apx C
import matplotlib.pyplot as plt
import numpy as np
import colorednoise as cn
from wetools import norm
import pyentrp.entropy as ent
import scipy.integrate as spi
from obspy.core.trace import Trace

r = 0.03                 # threshold
beta = 0                 # white noise
samples = int(1000/0.15) # number of samples to generate - means we have dt of 0.005

fname_start_no = 0

# Create original sinusoid
x = np.linspace(0, 1000, samples)
sinusoid = norm(0.5 * np.sin(2 * np.pi * x) + np.sin(3 * np.pi * x) + 0.7 * np.sin(3.6 * np.pi * x))

# Create different coefficient values
num_coeffs = 100
max_C = 300
C = np.linspace(0, max_C, num_coeffs)

# Do 100 times
for k in range(100):

    vals = [[], [], [], []]

    # Create noise
    n1 = norm(cn.powerlaw_psd_gaussian(beta, samples))
    N = np.random.normal(loc=0, scale=1, size=samples)

    # Set different bandpasses to test
    clrs = ['r', 'b', 'g', 'orange']
    bplow   = [0.01,  0.01, 1,   1.5]
    bpshigh = [2.0 ,  1.0 , 2.0, 2.0]

    # Loop through coefficient that is varying
    for coeff in C:
        # Using coeff = 50 for anything being held constant

        # Precursor effect
        front_width =  coeff/15
        taper = np.zeros(samples)
        taper[np.logical_and(x >= 5, x < 5+front_width)] = (1 - np.cos((x[np.logical_and(x >= 5, x < 5+front_width)] - 5) *  np.pi / (front_width))) / 2
        taper[np.logical_and(x >= 5+front_width, x < 6+ 2*front_width)] = 1

        # Coda Q effect:
        taper[x >= 6+ 2*front_width] = norm(np.exp(-x[x >= 6+ 2*front_width] /(1+ coeff))  )
        #taper[x >= 6+ 2*front_width] = norm(np.exp(-x[x >= 6+ 2*front_width] /(1+ 50))  )


        # Noise in taper
        n2 = norm(np.power(np.abs(N), 5 - 5*coeff/max_C  ))
        taper[x>=5] =  norm(  norm(taper[x>=5] *n2[x>=5])*taper[x>=5] )

        # Loop different bandpasses:
        for i in range(len(bplow)):

            # Noise in trace:
            timeseries =  norm((sinusoid + 2*n1*coeff/max_C)*taper)
            tr = Trace()
            tr.data = timeseries
            tr.stats.delta = x[2]-x[1]

            # Filter:
            tr = tr.filter(type='bandpass',freqmin=bplow[i], freqmax=bpshigh[i])

            timeseries = norm(tr.data)

            # Compute MSE and store
            MSE = ent.multiscale_entropy(time_series=timeseries[x>5], sample_length=2, tolerance=r, maxscale=200)
            m = np.arange(1, 201)
            vals[i].append(spi.trapz(MSE, m))


    # Convert to numpy array - rows = different coeff values, cols = different BPs
    # Save data:
    v = np.array(vals).T
    np.savetxt(fname=f"./synthetic_MSE_apx/data/combined_{k+fname_start_no}", X=v)
    print(f'Finished {k+fname_start_no}')



