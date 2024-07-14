from qc_funcs import *
from load_avg_Qc import calc_avg_QC

# Plot example Qc calculation (Appendix Figure A5)

# Plotting offset
offset = 5
CUTOFF = 150

# Choose parameter space simulation
dv = 'p-0.2'
a_list = ["1.5"]
l_list = [3]
stns   = "q"
path   = "../../data/processed/"

for l in l_list:
     for a in a_list:
        for stn in stns:

            # Load time data
            time = np.loadtxt(f"{path}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/time_data")
            dt = time[1]-time[0]

            # Load traces
            R = np.loadtxt(f"{path}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/M_0.ST{stn}.R")
            T = np.loadtxt(f"{path}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/M_0.ST{stn}.T")
            Z = np.loadtxt(f"{path}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/M_0.ST{stn}.Z")

            # Compute Qc and plot
            fig, ax, Qc_arr, time, energy = get_QC(time=time,
                                                    R=R,
                                                    T=T,
                                                    Z=Z,
                                                    convert_to_vel=True)


            # Array index for max energy (normalised)
            ind = np.where(energy == 1)[0][0]
            mean_orig = np.mean(Qc_arr[:, 2])
            std_orig = np.std(Qc_arr[:, 2])


            # Plot original decay:
            ax[2] = plot_decay_curve(ax=ax[2], mean=mean_orig, std=std_orig, time=time, ind=ind, freq=1, offset=offset)

            mean_reduced, std_reduced = calc_avg_QC(d=Qc_arr, upper_cutoff=CUTOFF)
            ax[3] = plot_decay_curve(ax=ax[3], mean=mean_reduced, std=std_reduced, time=time, ind=ind, freq=1,
                                     offset=offset)

            ax[2].set_title(f"{mean_orig} +- {std_orig}: offset: {offset}")
            ax[3].set_title(f"{mean_reduced} +- {std_reduced}: offset: {offset}")

            fig.suptitle(f"a = {a}  l = {l}  stn = {stn}")
            for i in range(4):
                ax[i].set_ylabel("Norm. Energy")

            ax[-1].set_xlabel("Time [s]")

            print(mean_orig, std_orig)
            print(mean_reduced, std_reduced)

plt.savefig(fname=f"example_Qc_plot.pdf", dpi='figure', format='pdf')
