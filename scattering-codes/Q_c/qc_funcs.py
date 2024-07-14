# Some functions used for Qc calculations
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as ss
from wetools import norm

def get_Qc_indices(q, energy):
    # Seeks time series index with max energy
    peak_ind = np.where(energy ==0)
    if len(peak_ind) != 1:
        print(f"Warning: More than one peak index (log norm energy = 0). Actual number: {len(peak_ind)}")
    last_ind = np.where(energy > q)[-1][-1]
    peak_ind = peak_ind[0][0]
    return [peak_ind, last_ind]


def moving_avg(f, time=[], half_window=1000, convert_t=False):
    # Computes moving average time series based on the half window length
    # If convert_t = True will also generate corresponding time array
    # Initialise mov_avg array:
    avg = np.zeros(int(len(f) - 2 * half_window))

    # Moving-average calculations cuts off the first and last number of elements equal to the value of 'half_window'
    for ts_index in range(half_window, int(len(f) - half_window)):
        avg[ts_index - half_window] = 1 / (half_window * 2) * np.sum(f[ts_index - half_window:ts_index + half_window])

    # If convering a time series eg u(t), may wish to slice time array to same size:
    if convert_t == True:
        try:
            t_avg = time[half_window:int(len(f) - half_window)]
            return t_avg, avg
        except time == []:
            print('Error: No time array inputted')
    else:
        return avg


def calc_grad(u, dt):
    # 4th order first derivative of time: v = \partial_t u = -(1/12)U[t+2] + (2/3)U[t+1] - (2/3)U[t-1] + (1/12)U[t-2]
    a1 = 1/12.0
    a2 = 2.0/3.0
    B = ss.diags([a1, -a2, a2, -a1], offsets=[-2, -1, 1, 2], shape=(len(u), len(u)) )
    v = (B*u)/dt
    v = v[2:-2]
    return v

def get_QC(time, R, T, Z, convert_to_vel, n=0):
    #Computes Qc based on displacement/velocity data
    # and an 'n' value (1/t^n)
    dt = time[1]-time[0]

    # Ensure traces are in velocity
    if convert_to_vel:
        time = time[2:-2]
        v_R = calc_grad(R,dt)
        v_T = calc_grad(T,dt)
        v_Z = calc_grad(Z,dt)
    else:
        print("Warning: not converting to velocity")
        time = time
        v_R = R
        v_T = T
        v_Z = Z

    # Sum energy from all orthogonal directions
    RTZ =  [v_R,v_T,v_Z]
    energy = RTZ[0]*0
    for c in RTZ:
        energy += c**2


    fig_result, ax_result = plt.subplots(5, sharex=True, figsize=(10,7))
    # plot norm vertical velocity
    ax_result[4].plot(time, norm(v_Z), 'k')

    #fig, ax = plt.subplots(4, sharex=True)
    #ax[0].set_xlim([0, 50])
   # for i in range(3):
        #ax[i].plot(time, RTZ[i], 'k')

    # Loop through various q cutoff values
    q_cutoffs = np.arange(-2.5, -1.0, 0.2)
    freq = 1
    data = []

    # plot energy data on 4 subplots for other stuff to be overlayed
    for i in range(4):
        ax_result[i].plot(time, norm(energy), 'k')
        ax_result[i].set_ylim([0, 1])

    # Normalise the energy and multiply by t^n
    energy = norm(energy * (time**n))

    # Loop through window size:
    for window in range(1, 101): #101
        # Plotting energy:
        t_avg, energy_avg = moving_avg(energy, time=time, half_window=window, convert_t=True)



        # Isolate the region that falls within log(energy) = 0 and
        # log(enegy) = q where q is 4 or 4.5 or 5 - undecided - may avg. all
        # Find index for which log energy average is 0 (ie peak)
        log_energy = np.log(norm(energy_avg))

        for q in q_cutoffs:
            ind = get_Qc_indices(q, log_energy)
            # Plot smoothed log energy
            #ax[4].plot(t_avg[ind[0]:ind[1] + 1], log_energy[ind[0]:ind[1] + 1])

            # Plot section sampled for Qc on seismograms
            #for i in range(3):
               #ax[i].plot(time[ind[0]:ind[1] + 1], RTZ[i][ind[0]:ind[1] + 1])


            # Calculate Qc:
            qc_trace_time = t_avg[ind[0]:ind[1] + 1]
            qc_trace      = log_energy[ind[0]:ind[1] + 1]
            p = np.polyfit(qc_trace_time, qc_trace, deg=1)
            Qc = (-2*np.pi*freq)/p[0]

            # Plot decay on original energy vs time
            t_exp = t_avg[ind[0]:]
            ax_result[1].plot(t_exp, norm(np.exp(-(2*np.pi*freq*t_exp)/Qc)), color='purple', alpha=0.05 )

            data.append([window, q, Qc])


    D = np.array(data)
    return fig_result, ax_result, D, time, norm(energy)











def get_QC_interpolate(time, R, T, Z, convert_to_vel, glitches, starttime):

    dt = time[1]-time[0]

    # Ensure traces are in velocity
    if convert_to_vel:
        time = time[2:-2]
        v_R = calc_grad(R,dt)
        v_T = calc_grad(T,dt)
        v_Z = calc_grad(Z,dt)
    else:
        print("Warning: not converting to velocity")
        time = time
        v_R = R
        v_T = T
        v_Z = Z

    # FIX GLITCHES:
    for gl in range(len(glitches)):
        gstart = glitches[gl][0]
        gend = glitches[gl][1]

        mask = np.logical_and(time>gstart+starttime, time<gend+starttime)
        v_R[mask] = 0
        v_T[mask] = 0
        v_Z[mask] = 0

    RTZ =  [v_R,v_T,v_Z]





    # Sum energy from all orthogonal directions
    energy = RTZ[0]*0
    for c in RTZ:
        energy += c**2


    # Now interpolate glitched regions:
    for gl in range(len(glitches)):
        gstart = glitches[gl][0]
        gend = glitches[gl][1]

        # get time value index for the beginning and end of glitches
        i1 = np.where(time==time[time <= gstart+starttime][-1])[0][0]
        i2 = np.where(time==time[time <= gend+starttime][-1])[0][0]

        # Now calculate line between values at i1 and i2:
        e1 = energy[i1]
        e2 = energy[i2]
        grad = (e2-e1)/(time[i2]-time[i1])



        energy[i1:i2 + 1] = grad* (time[i1:i2 + 1] - time[i1]) + e1

    fig_result, ax_result = plt.subplots(5, sharex=True, figsize=(10,7))

    ax_result[4].plot(time, norm(v_Z), 'k')

    #fig, ax = plt.subplots(4, sharex=True)
    #ax[0].set_xlim([0, 50])
   # for i in range(3):
        #ax[i].plot(time, RTZ[i], 'k')

    q_cutoffs = np.arange(-2.5, -1.0, 0.2)
    freq = 1


    data = []
    # Loop through window size:
    for i in range(4):
        ax_result[i].plot(time, norm(energy), 'k')
        ax_result[i].set_ylim([0, 1])

    energy = norm(energy)

    for window in range(1, 101): #101

        # Plotting energy:
        t_avg, energy_avg = moving_avg(energy, time=time, half_window=window, convert_t=True)

        #ax[3].plot(t_avg, norm(energy_avg), 'r')
        #ax[4].plot(time, np.log(norm(energy)), 'k')
        #ax[4].plot(t_avg, np.log(norm(energy_avg)), 'r')
        #ax[4].set_ylim([-7, 0])

        #plt.show()

        # Isolate the region that falls within log(energy) = 0 and log(enegy) = q where q is 4 or 4.5 or 5 - undecided - may avg. all
        # Find index for which log energy average is 0 (ie peak)
        log_energy = np.log(norm(energy_avg))

        for q in q_cutoffs:
            ind = get_Qc_indices(q, log_energy)
            # Plot smoothed log energy
            #ax[4].plot(t_avg[ind[0]:ind[1] + 1], log_energy[ind[0]:ind[1] + 1])

            # Plot section sampled for Qc on seismograms
            #for i in range(3):
               #ax[i].plot(time[ind[0]:ind[1] + 1], RTZ[i][ind[0]:ind[1] + 1])


            # Calculate Qc:
            qc_trace_time = t_avg[ind[0]:ind[1] + 1]
            qc_trace      = log_energy[ind[0]:ind[1] + 1]
            p = np.polyfit(qc_trace_time, qc_trace, deg=1)
            Qc = (-2*np.pi*freq)/p[0]

            # Plot decay on original energy vs time
            t_exp = t_avg[ind[0]:]
            ax_result[1].plot(t_exp, norm(np.exp(-(2*np.pi*freq*t_exp)/Qc)), color='purple', alpha=0.05 )

            data.append([window, q, Qc])


    D = np.array(data)
    return fig_result, ax_result, D, time, norm(energy)


def norm(x):
    return x/np.amax(np.abs(x))

def calc_exp_curve(t_exp, freq, Qc,n):
    return norm( (1/(t_exp**n))* np.exp(-(2 * np.pi * freq * t_exp) / Qc))

def plot_decay_curve(ax, mean, std, time, ind, freq=1, offset=0, n=0):
    exp = calc_exp_curve(t_exp=time[ind:], freq=freq, Qc=mean, n=n)
    ax.plot(time[ind:]-offset, exp, 'green')  # Plot mean
    exp_plus = calc_exp_curve(t_exp=time[ind:], freq=freq, Qc=mean - std, n=n)
    exp_minus = calc_exp_curve(t_exp=time[ind:], freq=freq, Qc=mean + std, n=n)
    ax.fill_between(time[ind:]-offset, y1=exp_minus, y2=exp_plus, color='green', alpha=0.5)

    ax.plot(time[ind:]-offset, exp_plus, 'blue')  # Plot mean
    ax.plot(time[ind:]-offset, exp_minus, 'blue')  # Plot mean

    return ax