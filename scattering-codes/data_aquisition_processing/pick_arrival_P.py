import numpy as np
stn_no = {"b": 0, "c": 1, "d": 2, "e": 3, "f": 4, "g": 5, "h": 6, "i": 7, "j": 8, "k": 9, "l": 10, "m": 11,
          "n": 12, "o": 13, "p": 14, "q": 15, "r": 16, "s": 17, "t": 18, "u": 19, "v": 20, "w": 21, "x": 22, "y": 23}
stn_list = "bcdefghijklmnopqrstuvwxy"

# Import filtered GF:
chl = 'R'
angle = 'M_0'
delta_v = 'p0.2'
out_dir = f"../../data/picks/p_arrival/{delta_v}/"

a_list = ["0.5", "1", "1.5", "2", "2.5"]
l_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

for l in l_list:
    for a in a_list[:l]:
        data_dir = f"../../data/processed/{delta_v}/{delta_v}_2hz_{l}_mfp_{a}_rad"

        t = np.loadtxt(f"{data_dir}/time_data")
        global_pick = []

        ctr = 0
        # Loop for each station
        for stn in stn_list:

            data = np.loadtxt(fname=f"{data_dir}/M_0.ST{stn}.{chl}")
            data = data/np.amax(np.abs(data))

            # Ensure that pick isnt within 2s time period (detected bc convolution artifact rather than arrival)
            min_time = t[t >= 2]
            data = data[t >= 2]
            pick = min_time[np.abs(data) >= 0.001][0]
            global_pick.append(pick)

            # Double check nothing is wrong with decreasing pick time for later stations
            if ctr > 0:
                assert global_pick[ctr] - (global_pick[ctr]-1) > 0
            ctr += 1


        out_str = f"{out_dir}/{l}_{a}"
        np.savetxt(fname=out_str, X=np.array(global_pick), fmt="%f")
        print(f"Saved {out_str}")
