import matplotlib.pyplot as plt
import numpy as np
from wetools import moving_avg, norm
import scipy.sparse as ss
# Pick the max energy of AxiSEM3D synthetic data
# - this may not be the exact arrival of the S phase
# - but is what we use for slicing off the qc calculations (when there is no moving window average)
# ----------------------------------------------------------------------------------------------------------------------

def calc_grad(u, dt):
    # 4th order first derivative of time:
    # v = \partial_t u = -(1/12)U[t+2] + (2/3)U[t+1] - (2/3)U[t-1] + (1/12)U[t-2]
    a1 = 1/12.0
    a2 = 2.0/3.0
    B = ss.diags([a1, -a2, a2, -a1], offsets=[-2, -1, 1, 2], shape=(len(u), len(u)) )
    v = (B*u)/dt
    v = v[2:-2]
    return v
# ----------------------------------------------------------------------------------------------------------------------

def get_Qc_indices(q, energy):
    peak_ind = np.where(energy ==0)
    if len(peak_ind) != 1:
        print(f"Warning: More than one peak index (log norm energy = 0). Actual number: {len(peak_ind)}")
    last_ind = np.where(energy > q)[-1][-1]
    peak_ind = peak_ind[0][0]
    return [peak_ind, last_ind]

# ----------------------------------------------------------------------------------------------------------------------
# Load trace:
dv = 'p-0.2'                                    # Pos/neg perturbation dataset
path = "../../data/processed/"                  # Path to processed data
out_dir =  f"../../data/picks/maxE_arrival/{dv}"   # Output path

# Set a value (radius) and loop through l values (separations)
a = "2.5"
for l in range(6, 11):
    # Load time and set timestep
    time = np.loadtxt(f"{path}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/time_data")
    dt = time[1]-time[0]

    # Loop each stations, compute velocity and then energy
    picks = []
    for stn in "bcdefghijklmnopqrstuvwxy":
        R = np.loadtxt(f"{path}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/M_0.ST{stn}.R")
        T = np.loadtxt(f"{path}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/M_0.ST{stn}.T")
        Z = np.loadtxt(f"{path}/{dv}/{dv}_2hz_{l}_mfp_{a}_rad/M_0.ST{stn}.Z")
        time = time[2:-2]
        v_R = calc_grad(R,dt)
        v_T = calc_grad(T,dt)
        v_Z = calc_grad(Z,dt)

        RTZ =  [v_R,v_T,v_Z]

        # Calculate energy (normalised)
        energy = RTZ[0]*0
        for c in RTZ:
            energy += c**2
        energy = norm(energy)

        ind_max = np.where(norm(energy)==1)[0][0]
        picks.append(time[ind_max])

    # Save picks to file
    out_path = f"{out_dir}/{l}_{a}"
    np.savetxt(fname=out_path, X=np.array(picks))
    print(f"Saved {out_path}")