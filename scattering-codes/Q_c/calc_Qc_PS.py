import numpy as np
import matplotlib.pyplot as plt
from qc_funcs import calc_grad, moving_avg, norm, get_Qc_indices

# Load traces for simulation parameters in list below:
dv     = 'p-0.2'
a_list = ["1.5"]
l_list = [3]

stns    = "q"
path    = "../../Scattering/data/processed/"
outpath = "../../data/Qc/"

for l in l_list:
     for a in a_list:
        for stn in stns:

            # load time data
            time = np.loadtxt(f"{path}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/time_data")
            dt = time[1]-time[0]

            # load channel data and convert to velocity
            R = np.loadtxt(f"{path}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/M_0.ST{stn}.R")
            T = np.loadtxt(f"{path}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/M_0.ST{stn}.T")
            Z = np.loadtxt(f"{path}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/M_0.ST{stn}.Z")
            time = time[2:-2]
            v_R = calc_grad(R,dt)
            v_T = calc_grad(T,dt)
            v_Z = calc_grad(Z,dt)

            # Compute energy
            RTZ =  [v_R,v_T,v_Z]
            energy = RTZ[0]*0
            for c in RTZ:
                energy += c**2

            # Plot channel data
            fig, ax = plt.subplots(5, sharex=True)
            for i in range(3):
                ax[i].plot(time, RTZ[i], 'k')

            # Vary q cutoff values
            q_cutoffs = np.arange(-2.5, -1.0, 0.2)
            freq = 1


            data = []
            # Loop through window size:
            for window in range(1, 101):

                # Plotting energy:
                t_avg, energy_avg = moving_avg(energy, time=time, half_window=window, convert_t=True)

                ax[3].plot(t_avg, norm(energy_avg), 'r')
                ax[3].plot(time, norm(energy), 'k')
                ax[3].set_ylim([0, 1])
                ax[4].plot(time, np.log(norm(energy)), 'k')
                ax[4].plot(t_avg, np.log(norm(energy_avg)), 'r')
                ax[4].set_ylim([-7, 0])


                # Isolate the region that falls within log(energy) = 0 and log(enegy) = q where q is 4 or 4.5 or 5 - undecided - may avg. all
                # Find index for which log energy average is 0 (ie peak)
                log_energy = np.log(norm(energy_avg))

                for q in q_cutoffs:
                    ind = get_Qc_indices(q, log_energy)
                    # Plot smoothed log energy
                    ax[4].plot(t_avg[ind[0]:ind[1] + 1], log_energy[ind[0]:ind[1] + 1])

                    # Plot section sampled for Qc on seismograms
                    for i in range(3):
                       ax[i].plot(time[ind[0]:ind[1] + 1], RTZ[i][ind[0]:ind[1] + 1])


                    # Calculate Qc:
                    qc_trace_time = t_avg[ind[0]:ind[1] + 1]
                    qc_trace      = log_energy[ind[0]:ind[1] + 1]
                    p = np.polyfit(qc_trace_time, qc_trace, deg=1)
                    Qc = (-2*np.pi*freq)/p[0]

                    # Plot decay on original energy vs time
                    t_exp = t_avg[ind[0]:]
                    ax[3].plot(t_exp, norm(np.exp(-(2*np.pi*freq*t_exp)/Qc)) )

                    data.append([window, q, Qc])

            # Save data:
            D = np.array(data)
            out_str = f"./{dv}/{dv}_{l}_{a}_{stn}"
            #np.savetxt(fname=out_str, X=D)
            print(f"Saved {out_str}")

plt.show()