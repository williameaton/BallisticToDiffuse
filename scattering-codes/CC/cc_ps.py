import matplotlib.pyplot as plt
import numpy as np
from calc_CC import calculate_CC

# Define simulation parameters:
dv   = 'p0.2'                              # Perturbation
chl  = "Z"                                 # Channel
stns = "bcdefghijklmnopqrstuvwxy"          # Station list
T    = 10                                  # Window time
ddir = "../../data/MSE_Pslices_250/"       # Data directory
odir = "../../data/CC/T=10_PROPERCUTOFF"   # Output directory
fmin = 0.1                                 # Min frequency
fmax = 2.0                                 # Max frequency

# List of a and l values
a_list = ["0.5"]
l_list = [1]

# Loop over a and l values
for l in l_list:
    for a in a_list:
        ctr = 0

        # Initialise mean_CC arrays to store mean values for all stations
        mean_CC = []

        # Loop over stations
        for stn in stns:

            # Define file name for simulation station, load sliced data + get timestep
            file = f"{ddir}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/ST{stn}.{chl}"
            sliced = np.loadtxt(fname=file)
            dt = sliced[1,0] - sliced[0,0]

            # Compute frequency of sliced data
            freq, C, Cm, nowindows = calculate_CC(sliced[:,1], dt=dt, T=T)


            # Generate plot for each station
            fig, ax = plt.subplots(2)

            # Slice upper freq
            maskhigh = freq <= 2.0
            fhigh  = freq[maskhigh]

            CC = C[:len(fhigh), :len(fhigh)]
            ax[1].set_title(f"Station: {stn}")
            ax[1].set_xlabel("Freq")
            ax[1].set_ylabel("Freq")
            ax[1].imshow(CC)


            # Need to use lower slicing properly
            indl = np.where(freq >= fmin)[0][0]
            indh = np.where(freq > fmax)[0][0]

            freq = freq[indl:indh]
            CCC = C[indl:indh, indl:indh]
            cm = np.mean(C[np.triu_indices_from(CCC, k=1)])
            print("Mean CC value within frequency band: ", cm)
            mean_CC.append(cm)

        ## Define output path:
        ## Uncomment to save values
        path = f"{odir}/{dv}/T_{T}_{l}_{a}.{chl}"
        #np.savetxt(fname=path, X=mean_CC, fmt="%f")
        #print(f"Saved {path}")

## Uncomment to show figures:
#plt.show()
