# Compute MSE for parameter space simulation:
# Here we compute MSE for region after max energy arrival
import numpy as np
import os
import pyentrp.entropy as ent
import scipy.integrate as spi

# Directories:
datadir = "../../data/MSE_Sslices_200"               # Input dir of time series
outdir  = "../../data/entropy/past_max_S_energy"     # Output dir for MSE

# Seismogram parameters:
dv   = "p-0.2"               # perturbation: either p0.2 or p-0.2
chl  = "T"                   # station channel
r    = 0.01                  # Threshold
maxm = 200                   # max m value
fmax = 2.0                   # Hz

# Station string
stns = "bcdefghijklmnopqrstuvwxy"


# Set m (scale) and f (fmax/m) arrays
m = np.arange(1,maxm+1)
f = fmax/m

# For looping over a and l
a_list = ["0.5", "1", "1.5", "2", "2.5"]
l_list = [2, 3, 4, 5, 6, 7, 8, 9, 10]

# loop l
for l in l_list:
    # some values dont have all l values:
    if l == 2 or l == 4 or l ==5:
        p = l-1
    else:
        p = l

     # loop a values
    for a in a_list[:p]:

        # output directories for raw MSE and integrated MSE
        out_dir_MSE    = f"{outdir}/MSE/r{r}/{dv}/{l}_{a}"
        out_dir_intMSE = f"{outdir}/intMSE/r{r}/{dv}/"

        # create directories if needed
        for path in [out_dir_intMSE, out_dir_MSE]:
            if os.path.isdir(path)==False:
                os.mkdir(path)

        # initialise for integration
        area = []

        # loop stations
        for stn in stns:
            # Load data
            fname = f"{datadir}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/ST{stn}.{chl}"
            data = np.loadtxt(fname)
            data = data[:,1]

            # Note we dont re-normalise the data as it is normalised before slicing
            MSE = ent.multiscale_entropy(time_series=data, sample_length=2, tolerance=r, maxscale=200)

            # Save the MSE:
            out_name_MSE = f"{out_dir_MSE}/ST_{stn}.{chl}"
            np.savetxt(out_name_MSE, MSE)
            print(f"Saved {out_name_MSE}")

            # Calculate the area under m-curve and freq-curve:
            Tm = spi.trapz(MSE, m)
            Tf = spi.trapz(MSE[::-1], f[::-1])

            # add integrated Tm and Tf data
            area.append([Tm, Tf])

        # Output data
        out_intMSE = f"{out_dir_intMSE}/{l}_{a}_{chl}"
        np.savetxt(fname=out_intMSE, X=area)
        print(f"Saved {out_intMSE}")