# Reproduce Fig A13 (supplementary material)
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(2, 3, sharex=True, sharey=True)

F = [[0.1, 2], [0.1, 0.6], [0.6, 1.1], [1.1, 1.75], [1.75, 2.0]]

# Loop through the frequencies
for M in range(len(F)):
    # Bandpass
    fmin = F[M][0]
    fmax = F[M][1]

    m = np.arange(1,201)

    if M < 3:
        i=0
        j = M
    else:
        i=1
        j = M-3

    # Load lunar data
    a1971 = np.loadtxt(f'./lunar/REV_MSE_1971a_14_Z_{fmin}_{fmax}')
    b1971 = np.loadtxt(f'./lunar/REV_MSE_1971b_12_Z_{fmin}_{fmax}')
    d1972 = np.loadtxt(f'./lunar/REV_MSE_1972_12_Z_{fmin}_{fmax}')

    # Load Martian
    if fmax ==2:
        fmax = '2.0'
    S0173 = np.loadtxt(f'./mars/JREV_FINAL_MSE_S0173_Z_{fmin}_{fmax}')
    S0235 = np.loadtxt(f'./mars/JREV_FINAL_MSE_S0235_Z_{fmin}_{fmax}')

    # Plot
    ax[i,j].plot(m, S0173, 'k',      label='S0173')
    ax[i,j].plot(m, S0235, 'purple', label='S0235')
    ax[i,j].plot(m, a1971, 'r',      label='1971a')
    ax[i,j].plot(m, b1971, 'g',      label='1971b')
    ax[i,j].plot(m, d1972, 'b',      label='1972')
    ax[i,j].set_title(f"{fmin} - {fmax}")

    ax[i,j].legend()

ax[0,0].set_ylim([0,1.5])
ax[0,0].set_xlim([1,100])
plt.show()
