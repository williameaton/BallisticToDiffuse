import numpy as np

def moving_avg(f, time=[], half_window=1000, convert_t=False):
    # =====================================================================================================================================
    # DESCRIPTION:
    # Function produces a moving-average time series of user-inputted time-series with window-size 2*spacing+1
    # In the case that u = u(t), user may want accompanying t time-series to be calculated
    # INPUTS:
    #    f [1D array]             - number of seconds
    #    time (opt) [int]         - coarse-grain scale
    #    time (opt) [int, array]  - default set to 'No' indicates time array doesnt need to be calculated; otherwise 1D time array/time-series
    # OUTPUTS:
    #    cgu [1D array]           - coarse-grained time-series
    #    cg_time (opt) [1D array] - accompanying time array for time-series
    # =====================================================================================================================================


    # Initialise mov_avg array:
    avg = np.zeros(int(len(f) - 2 * half_window))

    # Moving-average calculations cuts off the first and last number of elements equal to the value of 'half_window'
    for ts_index in range(half_window, int(len(f) - half_window)):
        avg[ts_index - half_window]= 1 / (half_window * 2) * np.sum(f[ts_index - half_window:ts_index + half_window])

    # If convering a time series eg u(t), may wish to slice time array to same size:
    if convert_t==True:
        try:
            t_avg = time[half_window:int(len(f) - half_window)]
            return t_avg, avg
        except time == []:
            print('Error: No time array inputted')

    else:
        return avg



def slice_by_time(t, d, t_low, t_high):
    # Get indices of values to slice by:
    # Lower bound:
    lb = np.where(t>=t_low)
    lb = lb[0][0]
    # Upper bound:
    ub = np.where(t<=t_high)
    ub = ub[0][-1]
    return t[lb:ub+1], d[lb:ub+1]


def norm(x):
    return x/np.amax(np.abs(x))





def load_int_MSE(path, dv, channel):
    all_MSE = []
    pos = []
    a_list = ["0.5", "1", "1.5", "2", "2.5"]
    l_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    for l in l_list:
        if l == 2 or l == 4 or l ==5:
            p = l-1
        else:
            p = l
        for a in a_list[:p]:
            # Load integrated MSE data:
            intMSE = np.loadtxt(f"{path}/{dv}/{l}_{a}_{channel}")
            all_MSE.append(intMSE)
            pos.append([float(l),float(a)])

    return np.array(pos), np.array(all_MSE)


def load_avg_Qc(path, dv, norm_factor=120):

    if dv=="p0.2":
        pos_Qc = [["0.5", "1"], ["1", "2"],  ["1.5", "3"], ["2", "4"], ["2.5", "5"] ]
    elif dv=="p-0.2":
        pos_Qc = [["0.5", "1"], ["1", "3"],  ["1.5", "3"], ["2", "5"], ["2.5", "6"] ]
    else:
        raise ValueError("dv must be p-0.2 or p0.2")

    data_reduced = []
    data_original = []
    for i in pos_Qc:

        a = i[0]
        for l in range(int(i[1]), 11):

            original = []
            reduced = []

            if l == 5 and a == "1.5" and dv=="p0.2":
                stns_list = "bcdefghijklmnopsuvxy"
            elif l == 4 and a == "1" and dv == "p0.2":
                stns_list = "bcdefghijklmnoprstuvwxy"
            else:
                stns_list = "bcdefghijklmnopqrstuvwxy"

            stn_no = 0
            for stn in stns_list:
                d = np.loadtxt(f"{path}/{dv}_{l}_{a}_{stn}")
                d_max = np.amax(d[:, 2])



                d = d[(d[:, 2] > 0)]
                d = d[(d[:, 2] < 150)]

                #if len(d[:,2])==0:
                #   print(f"l = {l}   a = {a}    station: {stn}")

                # Data is now in format of [m, q, Qc ]
                mean = np.mean(d[:, 2])
                std_dev = np.std(d[:, 2])

                data_original.append([float(l), float(a), stn_no, mean, std_dev])

                # Remove outliers
                d = d[(d[:, 2] < (std_dev) + mean)]
                d = d[(d[:, 2] > -(std_dev) + mean)]
                mean = np.mean(d[:, 2])
                std_dev = np.std(d[:, 2])

                data_reduced.append([float(l), float(a), stn_no, mean, std_dev])
                stn_no += 1

            #data_original.append(original)
            #data_reduced.append(reduced)

    # Convert to np array:
    data_original = np.array(data_original)
    data_reduced = np.array(data_reduced)

    data_original[:, 2:] = data_original[:, 2:]   / norm_factor
    data_reduced[:, 2:]   = data_reduced[:, 2:]   / norm_factor

    return data_original, data_reduced


# Converts obspy trace into x and y for mpl plotting
def obspy_gen_mpl(tr):
    x = np.linspace(0, tr.stats.npts*tr.stats.delta,  tr.stats.npts)
    y = tr.data
    return x,y