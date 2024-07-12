# Processes the synthetic data for stated perturbation and l values:
# Raw data is stored in data/raw and outputted to data/processed
# Processing contains decimation and resampling to 200 Hz, as well as bandpassing.


import numpy as np
import os
import obspy.core.trace as TR

# Perturbation can be positive 20 % "p0.2" or negative "p-0.2"
dv = "p0.2"
# Set l value (separation
l = 10
lst = [f"{dv}_2hz_{l}_mfp_0.5_rad", f"{dv}_2hz_{l}_mfp_1_rad", f"{dv}_2hz_{l}_mfp_1.5_rad",
       f"{dv}_2hz_{l}_mfp_2_rad", f"{dv}_2hz_{l}_mfp_2.5_rad"]

# Output directory
out_dir = f"../../data/processed/{dv}/"

# Decimation and resampling frequency
dec = 2
resample_freq = 200.0
chls = "RTZ"

# Loop through each 'a' value for this 'l'
for folder in lst[:l]:
    if dv in folder:
        master_dir = f"../../data/raw/{dv}/{folder}/"
        print(f"master dir: {master_dir}")

        # Create directory in 'processed' if necessary
        if os.path.isdir(f"{out_dir}/{folder}/")==False:
            os.mkdir(f"{out_dir}/{folder}/")

        # Load time data
        time = np.loadtxt(fname=f"{master_dir}/data_time.ascii")
        time = time[::dec]
        dt = time[1] - time[0]


        for file in os.listdir(master_dir):
            if "" in file and "M_0" in file:
                print(f"Opening {file}")
                data = np.loadtxt(fname=f"{master_dir}/{file}")
                for ch in range(3):
                    # Create trace with decimated data
                    tr = TR.Trace()
                    tr.stats.delta = dt
                    d = data[:,ch]
                    tr.data = d[::dec]

                    # Downsample data to constant sampling rate of 200 Hz:
                    tr.resample(resample_freq)

                    # Filter data:
                    tr.filter(type="bandpass", freqmin=0.01, freqmax=2.0)

                    out_str = f"{out_dir}/{folder}/{file[:-6]}.{chls[ch]}"
                    np.savetxt(fname=out_str, X=np.transpose([tr.data]))
                    print(f"Processed {out_str}")


    # Write time
    time = np.linspace(0, tr.stats.delta * (tr.stats.npts - 1), tr.stats.npts)
    time_out_str = f"{out_dir}/{folder}/time_data"
    np.savetxt(fname=time_out_str, X=np.transpose([time]), fmt="%f")

