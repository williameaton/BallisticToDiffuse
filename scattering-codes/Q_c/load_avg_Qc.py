import numpy as np

# Computes average Qc
def calc_avg_QC(d, upper_cutoff):

    values = d[:,2]
    d = values[np.logical_not(np.isnan(values))]
    #d_max = np.amax(d[:, 2])
    d = d[(d > 0)]
    d = d[(d < upper_cutoff)]

    # Data is now in format of [m, q, Qc ]
    mean = np.mean(d)
    std_dev = np.std(d)


    # Remove outliers
    d = d[(d < (std_dev) + mean)]
    d = d[(d > -(std_dev) + mean)]
    mean = np.mean(d)
    std_dev = np.std(d)

    return mean, std_dev



# Loads average Qc
def load_avg_Qc(path, dv, upper_cutoff, norm_factor=1):

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
                stns_list = "bcdefghijklmnopsuvx"
            elif l == 4 and a == "1" and dv == "p0.2":
                stns_list = "bcdefghijklmnoprstuvwx"
            else:
                stns_list = "bcdefghijklmnopqrstuvwx"

            stn_no = 0
            for stn in stns_list:
                d = np.loadtxt(f"{path}/{dv}_{l}_{a}_{stn}")

                mean, std_dev = calc_avg_QC(d, upper_cutoff=upper_cutoff)
                data_original.append([float(l), float(a), stn_no, mean, std_dev])
                data_reduced.append([float(l), float(a), stn_no, mean, std_dev])
                stn_no += 1

    # Convert to np array:
    data_original = np.array(data_original)
    data_reduced = np.array(data_reduced)

    data_original[:, 2:] = data_original[:, 2:]   / norm_factor
    data_reduced[:, 2:]   = data_reduced[:, 2:]   / norm_factor


    return data_original, data_reduced



