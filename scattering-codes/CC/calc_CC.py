import matplotlib.pyplot as plt
from windows import create_windows
import numpy as np

def calculate_CC(data, dt, T, window_overlap=0):
    # Calculates correlation coefficient matrix for a 1D time-series
    # Inputs:
    #   data [array-like or numpy] :    1D time series
    #   dt   [float]               :    timestep of data
    #   T    [float]               :    window size in seconds
    #   window_overlap [float]     :    percentage overlap of sample windows
    # Returns:
    #   freq_out [array-like]      :    Array of frequencies for which CC is calculated
    #   C                          :    Correlation coefficient matrix
    #   C_mean                     :    Mean of lower triangluar of matrix (probably ignore this)
    #   no_windows                 :    number of windows generated


    # create windows of noise:
    windows, no_windows = create_windows(data=data, T=T, dt=dt, O=window_overlap, padded=True, summary=True)

    # Get length of resultant power array and initialise phi array:
    freq        = np.fft.fftfreq(len(windows[0,:]), dt)
    freq_out    = freq[freq > 0]
    phi_len     = len(freq[freq > 0])

    phi = np.zeros((no_windows, phi_len), dtype=np.cdouble)

    for i in range(no_windows):
        power = np.fft.fft(windows[i,:])
        phi[i,:] = power[freq > 0]


    # Create phi^2:
    phi_sq   = np.square(np.abs(phi))
    E_phi_sq = np.array([np.mean(phi_sq, axis=0)])
    E_phi_phi = np.zeros((phi_len, phi_len), dtype=np.cdouble)


    for i in range(no_windows):
        # For each window calculate the product of phi and phi * --> creates matrix
        # Need an element-wise average of these matrices so add them up and
        phi_i = np.array([phi[i,:]])
        E_phi_phi += np.multiply(np.transpose(phi_i), np.conjugate(phi_i))


    E_phi_phi = np.square(np.abs(E_phi_phi/no_windows))

    C = np.divide(E_phi_phi, np.multiply(np.transpose(E_phi_sq), E_phi_sq))

    C_mean = np.mean(C[np.triu_indices_from(C, k=1)])

    return freq_out, C, C_mean, no_windows





if __name__ == "__main__":
    import colorednoise as cn
    # Example using colorednoise packages from:
    # https://github.com/felixpatzelt/colorednoise/blob/master/colorednoise.py
    N = 70006
    T = 10
    noise = cn.powerlaw_psd_gaussian(0, N) # White noise (beta = 0)

    f, C, CCm, windows = calculate_CC(data=noise, dt=0.005, T=T)

    fig, ax = plt.subplots()
    ax.imshow(C, origin='lower', extent=(min(f),  max(f), min(f), max(f)))
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Frequency [Hz]')
    plt.show()