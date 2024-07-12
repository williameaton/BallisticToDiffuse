# Script reads in processed synthetic data (from data/processed/...) and slices it.
# It slices the time-series starting at the maximum energy pick value (data/picks/maxE_arrival)
# and ending 200 seconds later. The data is outputted to data/MSE_Sslices_200.
import os
import numpy as np

def norm(x):
    return x/np.amax(np.abs(x))

def slice_by_time(t, d, t_low, t_high):
    # Get indices of values to slice by:
    # Lower bound:
    lb = np.where(t>=t_low)
    lb = lb[0][0]
    # Upper bound:
    ub = np.where(t<=t_high)
    ub = ub[0][-1]
    return t[lb:ub+1], d[lb:ub+1]


# I cant believe I used to code like this:
stn_list = "bcdefghijklmnopqrstuvwxy"
stn_no = {"b": 0, "c": 1, "d": 2, "e": 3, "f": 4, "g": 5, "h": 6, "i": 7, "j": 8, "k": 9, "l": 10, "m": 11,
          "n": 12, "o": 13, "p": 14, "q": 15, "r": 16, "s": 17, "t": 18, "u": 19, "v": 20, "w": 21, "x": 22, "y": 23}

# Synthetic dataset - perturbation either positive (p0.2) or negative (p-0.2)
dv = "p-0.2"
cutoff = 200  # Seconds

# List of a values
a_list = ["0.5", "1", "1.5", "2", "2.5"]
# List of l values
l_list = [5]

for l in l_list:
    for a in a_list[:l-1]:
        # Path strings
        data_dir = "../../data/"
        sim = f"{dv}_2hz_{l}_mfp_{a}_rad"

        # Output path:
        out_path = f"{data_dir}/MSE_Sslices_{str(cutoff)}/{dv}/{sim}"

        # Create output folder if it doesnt exist
        if os.path.isdir(out_path) == False:
            os.mkdir(out_path)

        # Load pick of maximum energy arrival
        pick = np.loadtxt(f"{data_dir}/picks/maxE_arrival/{dv}/{l}_{a}")

        # Loop all stations
        for stn in stn_list:
            # Loop all channels
            for chl in "RTZ":
                # Load data
                t = np.loadtxt(f"{data_dir}/processed/{dv}/{sim}/time_data")
                d = np.loadtxt(f"{data_dir}/processed/{dv}/{sim}/M_0.ST{stn}.{chl}")
                d = norm(d)

                # Slice time series from pick value to pick value + 200 seconds
                ts, ds = slice_by_time(t=t, d=d, t_low=pick[stn_no[stn]], t_high=pick[stn_no[stn]]+cutoff)

                np.savetxt(fname=f"{out_path}/ST{stn}.{chl}", X=np.transpose(np.array([ts, ds])), fmt="%f" )
                print(f"Saved: {out_path}/ST{stn}.{chl}")
