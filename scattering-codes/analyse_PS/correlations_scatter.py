# Plots scattering correlations for Supplementary figures 8 - 11

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.append('../../plotting/')
from animate_PS import load_corr_data
sys.path.append('../Q_c/');
from load_avg_Qc import load_avg_Qc

# --------------------------------- SOME FUNCTIONS  ---------------------------------


def calc_vol_sa_etc(A,L,stn):
    # Computes volume and surface areas based on formula  (a/l)**3 , (a^2/l^3) etc
    # src at 20 km depth, 5 km station separation
    distance = ((20)**2 +   ((1+ stn)*5)**2)**0.5
    vol = ((A / L) ** 3)/0.125                      # Volume effect
    sa = (A ** 2 / L ** 3)/0.25                     # Shape effect
    al = L/(4-A)/6.666666666666667                  # Morphology effect

    return vol, sa, distance, al

# Function to load integrated MSE
def load_int_MSE(path, dv, channel='R', simlist=None):
    all_MSE = []
    pos = []

    if simlist==None:
        a_list = ["0.5", "1", "1.5", "2", "2.5"]
        l_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        for l in l_list:
            l = int(l)
            if l == 2 or l == 4 or l ==5:
                p = l-1
            else:
                p = l
            for a in a_list[:p]:
                # Load integrated MSE data:
                intMSE = np.loadtxt(f"{path}/{dv}/{l}_{a}_{channel}")
                all_MSE.append(intMSE)
                pos.append([float(l),float(a)])
    else:
        s = np.array(simlist)
        print("Loading MSE with SIMLIST ")
        a_list = s[:,0]
        l_list = s[:,1]

        for i in range(len(a_list)):
            a = a_list[i]
            l = l_list[i]

            # Load integrated MSE data:
            intMSE = np.loadtxt(f"{path}/{dv}/{l}_{a}_{channel}")
            all_MSE.append(intMSE)
            pos.append([float(l), float(a)])

    return np.array(pos), np.array(all_MSE)

# ------------------------------ SET SOME PLOTTING PARAMETERS ------------------------------

outname = 'correlations_scatter.pdf'

# Create figures
fig_a, ax_a       = plt.subplots(2,3, figsize=(7,5), sharex=True)
fig_l, ax_l       = plt.subplots(2,3, figsize=(7,5), sharex=True)
fig_vol, ax_vol   = plt.subplots(2,3, figsize=(7,5), sharex=True)
fig_surf, ax_surf = plt.subplots(2,3, figsize=(7,5), sharex=True)

for f in [fig_surf, fig_l, fig_a, fig_vol]:
    f.set_tight_layout(True)

# Colours for each parameter - MSE, Qc, CC
clrs = ['purple', 'teal', 'coral']

xlims_for_lbf = [0, 80]

# Loop through positive and negative data
iperturb = 0
for dv in['p0.2', 'p-0.2']:

    stn_str = "bcdefghijklmnopqrstuvwx"
    no_stns = len(stn_str)

    # Initialise empty stations
    dataCC  = [];  dataCC_LA  = []
    dataQC  = [];  dataQC_LA  = []
    dataMSE = [];  dataMSE_LA = []

    # Load data
    mse_pos, intMSE = load_int_MSE(path=f"../../data/entropy/intMSE/r0.01/", dv=dv)
    pos_CC, CC      = load_corr_data(path=f"../../data/CC/T=10_PROPERCUTOFF/{dv}",
                                     channel="Z", no_stations=no_stns, T=10)
    orig, red       = load_avg_Qc(path=f"../../data/Qc/{dv}/" , upper_cutoff=150, dv=dv, norm_factor=1)

    # Select stations for MSE:
    intMSE = intMSE[:,:no_stns, :]

    # format Qc:
    print(np.shape(red)[0])
    for i in range(np.shape(red)[0]):
        L = float(red[i, 0])
        A = float(red[i, 1])
        stn = float(red[i, 2])
        vol, sa, dist, al  = calc_vol_sa_etc(A,L,stn)
        dataQC.append([vol, sa, al, dist, red[i, 3]])

        dataQC_LA.append([L, A, dist, red[i, 3]])

    # format CC:
    print(np.shape(pos_CC)[0])
    for i in range(np.shape(pos_CC)[0]):
        L = float(pos_CC[i, 0])
        A = float(pos_CC[i, 1])
        for stn in range(no_stns):
            vol, sa, dist, al  = calc_vol_sa_etc(A,L,stn)

            dataCC.append([vol, sa, al, dist, CC[i, stn]])

            dataCC_LA.append([L,A, dist, CC[i, stn]])

    # format intMSE:
    print(np.shape(intMSE)[0])
    for i in range(np.shape(intMSE)[0]):
        L = float(mse_pos[i, 0])
        A = float(mse_pos[i, 1])
        for stn in range(no_stns):
            vol, sa, dist, al  = calc_vol_sa_etc(A,L,stn)
            # EDIT HERE TO FLIP BETWEEN INT M AND INT F - 0 for INTM 1 for INTF
            dataMSE.append([vol, sa, al, dist, intMSE[i, stn, 0]])
            dataMSE_LA.append([L, A, dist, intMSE[i, stn, 0]])

    # Convert to np arrays
    dataMSE    = np.array(dataMSE)
    dataCC     = np.array(dataCC)
    dataQC     = np.array(dataQC)

    dataMSE_LA = np.array(dataMSE_LA)
    dataCC_LA  = np.array(dataCC_LA)
    dataQC_LA  = np.array(dataQC_LA)

    # Set data arrays for looping
    data   = [dataMSE, dataQC, dataCC]          # Holds the volume and surface area
    dataLA = [dataMSE_LA, dataQC_LA, dataCC_LA] # Holds the a and l
    print('Plotting as numbers of wavelengths. Max is 116 km / (lambda = 1.5 km) ')


    # --------------- DATA IS NOW LOADED AND FORMATTED TIME TO PLOT ---------------
    # Plot data points
    for y in range(3):
        # Loop through CC, Qc, MSE
        for x in range(2):
            # If x = 0 then volume          and   l
            # If x = 1 then surface area    and   a

            # Determine axis, labels etc
            if x==0:
                axsolo = ax_l
                axcomb = ax_vol
                soloxlab = r'$\frac{l}{\lambda} \times \frac{L}{\lambda}$'
                combxlab = r'$\frac{a^3}{l^3} \times \frac{L}{\lambda}$'
            else:
                axsolo = ax_a
                axcomb = ax_surf
                soloxlab = r'$\frac{a}{\lambda} \times \frac{L}{\lambda}$'
                combxlab = r'$\frac{a^2}{l^3} \times \frac{L}{\lambda}$'


            # --------- Plot for Surface area or Volume ---------
            # Loop through [Volume, SA] * Station on x axis
            # x = 0: vol, x=1: sa
            # 1.5 km is wavelength
            X = data[y][:,x] * (data[y][:,-2])/1.5 # x Divide by wavelength in km
            Y = data[y][:,-1]

            # Calculate line of best fit polynomial
            p = np.polyfit(X, Y, deg=1)

            # Plot points:
            axcomb[iperturb, y].plot(X, Y, 'o', color=clrs[y], alpha=0.3)

            # Compute line of best fit for plot
            x_lbf = np.array(xlims_for_lbf)
            y_lbf = p[1] + p[0]*x_lbf

            # Plot lbf:
            axcomb[iperturb, y].plot(x_lbf, y_lbf, 'k', alpha=0.5)
            axcomb[iperturb, y].set_title(np.around(np.corrcoef(X, y=Y)[0,1],4 ))

            # Add x label
            axcomb[1, y].set_xlabel(combxlab)


            # -------- Plot of L or A against values:  -------
            # Compute xaxis (l or a times station distance)
            # stn distance divided by wavelength 1.5 km/s
            x_lora = dataLA[y][:,x] *  ((data[y][:,-2])/1.5 )
            D = data[y][:,-1]

            # plot
            axsolo[iperturb,y].plot(x_lora, D, 'o', color=clrs[y], alpha=0.3)

            # Add correlation values
            axsolo[iperturb,y].set_title(np.around(np.corrcoef(x_lora, y=D)[0, 1], 4))

            # Add label
            axsolo[1, y].set_xlabel(soloxlab)

    iperturb+=1

# Style figure in a VERY manual way:

# Set y limits for l and a plots:
ax_l[0, 0].set_xlim([0, 800])
ax_a[0, 0].set_xlim([0, 200])
ax_vol[0, 0].set_xlim([0, 80])
ax_surf[0, 0].set_xlim([0, 80])

for axi in [ax_l, ax_a, ax_surf, ax_vol]:
    for j in range(2):
        axi[j,0].set_ylim([0,310])
        axi[j,1].set_ylim([0,115])
        axi[j,2].set_ylim([1,0])


# save figures to a single file:
p = PdfPages(outname)
for fig in [fig_a, fig_l, fig_vol, fig_surf]:
    # and saving the files
    fig.savefig(p, format='pdf')

    # close the object
p.close()