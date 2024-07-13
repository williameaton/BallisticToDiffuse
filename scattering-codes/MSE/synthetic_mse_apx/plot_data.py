# Generates Apx C plots of synthetic time series

import numpy as np
import matplotlib.pyplot as plt

typ = 'combined_blue'                        # Dataset - see options in ./data subdirectory
name = f'correct_threshold/{typ}/{typ}_'     # # Set file name:
n_samples = 100                              # How many samples
n_coeffs  = 100                              # How many coefficients
clrs = ['r', 'g', 'b', 'orange']             # Plot colours

# Plotting params
min_mse = 0
max_mse = [10, 25, 3, 2]

# Create figure
fig, ax = plt.subplots(4,2, figsize=(12, 5))
fig.suptitle(f"{name} holding everything else at s=50")

# Loop through bandpasses
for bp in range(4):
    ex = []

    # Loop for each sample
    for k in range(n_samples):
        # Load
        data = np.loadtxt(f'./data/{name}{k}')

        # Plot individual data
        ax[bp, 0].plot( np.arange(n_coeffs)+1 , data[:,bp], color=clrs[bp], alpha=0.25)
        ax[bp, 0].set_xlim([1, n_coeffs])
        # Store for calculating mean/std etc etc
        ex.append(data[:,bp])


    # Convert to np array and initialise G for storing gaussians
    ex = np.array(ex)
    G = []

    # For each coefficient (s value sampled) get the mean/std of all the n_samples
    for i in range(n_coeffs):
        n_gauss = 1000
        mean = np.mean(ex[:,i])
        std  = np.std(ex[:,i])

        # Plot gaussian on line version
        x = np.linspace(min_mse, 50, n_gauss)
        g = 1/(std * (2*np.pi )**0.5 ) * np.exp(- (x-mean)**2/(2*(std**2))  )
        #ax[0].plot(g*0.2 +i + 1, x, 'k-', alpha=0.5)

        # Store for imshow plot
        G.append(g[:])

    # Format for imshow plot:
    G = np.array(G).T
    ax[bp, 1].imshow(G, aspect='auto', origin='lower', interpolation='none', extent=(1,n_coeffs, min_mse, max_mse[bp]),cmap='magma_r')
    ax[bp, 1].set_xlabel('No of coeff sample (from s=0 to s=300')
    ax[bp, 1].set_ylabel('Integrated MSE')

    ax[bp, 0].set_ylim([min_mse, max_mse[bp]])

# Output plots:
plt.savefig(fname=f'./figures/mseapx_{typ}.pdf', dpi='figure', format='pdf')
