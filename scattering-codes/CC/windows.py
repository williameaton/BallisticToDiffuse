# Function to create windows from a time series:
import numpy as np

def create_windows(data, T, dt, O, padded=True, summary=False):
    # Breaks up 1D time series into windows
    # Inputs:
    #   data [array-like or numpy] :    1D time series
    #   T    [float]               :    window size in seconds
    #   dt   [float]               :    timestep of data
    #   O    [float]               :    percentage overlap of sample windows
    #   padded [boolean]           :    if True, each window is doubled in length for padding with zeros
    #   summary [boolean]          :    print summary of windows
    # E.g. the following create_windows(data=d, T=100, dt=0.25, O=0.1)
    # Would mean produce windows for timeseries d of length 100 where sample spacing is 0.25 (ie datapoints per window)
    # and seperate each window by 10 % of T


    # Calculate number of windows required
    N = len(data)
    el_per_window = int(np.floor(T/dt))
    el_overlap    = int(O*el_per_window)
    no_win = int(np.floor( (N + el_overlap)/(el_per_window + el_overlap) ))

    # Initialise array of windows:
    if padded == True:
        windows = np.zeros((no_win, el_per_window*2))
    elif padded == False:
        windows = np.zeros((no_win, el_per_window))
    else:
        raise ValueError("padded must be boolean")



    for i in range(no_win):
        windows[i, :el_per_window] = data[(el_per_window+el_overlap)*i : (el_per_window+el_overlap)*i + el_per_window]

    if summary:
        print(f"Length of data:      {N}")
        print(f"Elements per window: {el_per_window}")
        print(f"Elements in overlap: {el_overlap}")
        print(f"Number of windows  : {no_win}")
        print(f"Using padding?     : {padded}")
    return windows, no_win






if __name__ == "__main__":
    from create_noise import gen_noise
    import matplotlib.pyplot as plt

    n = 2000
    DT = 0.1
    noise = gen_noise(n, color="white")
    t     = np.arange(n)*DT


    time = create_windows(t, T=13, dt=DT, O=0.1, padded=True)
    windows = create_windows(noise, T=13, dt=DT, O=0.1, padded=True)

    fig, ax = plt.subplots()
    ax.plot(t, noise, 'k-')
    for i in range(len(time[:,0])):
        ax.plot(time[i,:], windows[i,:], 'r-')

    plt.show()