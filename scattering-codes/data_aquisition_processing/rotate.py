# Functions to rotate stream/traces
# Following the method used by Jeroen Tromp in NMSYNG (Dahlen and Tromp 98 eqn 10.3)
# Note that I use 'TPZ' instead of 'RTZ' in convention with the 'Theta, Phi, Z' of DT98.
# This is because 'radial' from a spherical harmonic point of view indicates the vertical direction where as Theta
# ('commonly 'radial') is src-receiver direction and Phi is in the transverse direction.

import numpy as np

def rotate_stream(stream,  method, src, stn, geoco=1, overwrite_channel_names=False, invert_p=True, invert_t=True):
    # Get station coordinates:
    lat_stn = stn[0]
    lon_stn = stn[1]

    pi = np.pi
    lat_src = src[0]
    lon_src = src[1]

    # Convert geographic decimal --> geocentric radian values (if geoco==1 then geographic == geocentric):
    # Note here that theta is the CO-latitude
    theta_src = (pi/2) - np.arctan(geoco * np.tan( lat_src * pi/180 ))  # Source colatitude
    theta_stn = (pi/2) - np.arctan(geoco * np.tan( lat_stn * pi/180 ))  # Station colatitude

    phi_src   = lon_src*pi/180                                          # Source longitude
    phi_stn   = lon_stn*pi/180                                          # Station longitude

    # Calculate epicentral distance $\Theta$ (scalar):
    dist = np.arccos(np.cos(theta_stn)*np.cos(theta_src) + np.sin(theta_stn)*np.sin(theta_src)*np.cos(phi_stn - phi_src))

    rot1 = (1/np.sin(dist)) * \
           (np.sin(theta_stn)*np.cos(theta_src)  -  np.cos(theta_stn)*np.sin(theta_src)*np.cos(phi_stn - phi_src))

    rot2 = (1/np.sin(dist)) * (np.sin(theta_src)*np.sin(phi_stn - phi_src))

    # Conversion from RTP --> ZNE (where R=Z, 2D rotation) appears to use the following matrix:
    #   [N, E]' = [-rot1, -rot2; -rot2, rot1][T, P]' where T and P are theta, Phi
    #   Below we shall name the rotation matrix Q:
    # Hence to get the T and P matrix we should be multiplying [N,E] by the inverse of Q:
    Q    = np.array([[-rot1, -rot2], [rot2, -rot1]])
    Qinv = np.linalg.inv(Q)

    if method == "NE->TP":
        N = stream.select(component="N")[0].data
        E = stream.select(component="E")[0].data
        data_NE = np.array([N,E])
        data_TP = np.matmul(Qinv, data_NE)

        # Now writing back to stream:
        old_chls = ["N", "E"]
        new_chls = ["T", "P"]
        for i in range(2):
            if np.logical_and(new_chls[i] == "P", invert_p==True):
                data_TP[i, :] = data_TP[i,:]*(-1)
            if np.logical_and(new_chls[i] == "T", invert_t==True):
                data_TP[i, :] = data_TP[i,:]*(-1)

            stream.select(component=old_chls[i])[0].data = data_TP[i,:]
            if overwrite_channel_names:
                old_chl_name = stream.select(component=old_chls[i])[0].stats.channel

                if old_chl_name.find('E') != -1:
                    # Channel has E in it:
                    index = old_chl_name.find('E')
                    new_name = old_chl_name[:index] + 'P' + old_chl_name[index+1:]

                if old_chl_name.find('N') != -1:
                    # Channel has N in it:
                    index = old_chl_name.find('N')
                    new_name = old_chl_name[:index] + 'T' + old_chl_name[index+1:]

                stream.select(component=old_chls[i])[0].stats.channel = new_name

    else:
        raise ValueError("Currently method must be NE->TP")



# Following the method used by JT in NMSYNG (DT98 eqn 10.3)
def rotate_trace(traceN, traceE, method, src, stn, geoco=1):
    # Get station coordinates:
    lat_stn = stn[0]
    lon_stn = stn[1]

    pi = np.pi
    lat_src = src[0]
    lon_src = src[1]


    # Convert geographic decimal --> geocentric radian values (if geoco==1 then geographic == geocentric):
    # Note here that theta is the CO-latitude
    theta_src = (pi/2) - np.arctan(geoco * np.tan( lat_src * pi/180 ))  # Source colatitude
    theta_stn = (pi/2) - np.arctan(geoco * np.tan( lat_stn * pi/180 ))  # Station colatitude

    phi_src   = lon_src*pi/180                                               # Source longitude
    phi_stn   = lon_stn*pi/180                                               # Station longitude

    # Calculate epicentral distance $\Theta$ (scalar):
    dist = np.arccos(np.cos(theta_stn)*np.cos(theta_src) + np.sin(theta_stn)*np.sin(theta_src)*np.cos(phi_stn - phi_src))

    rot1 = (1/np.sin(dist)) * \
           (np.sin(theta_stn)*np.cos(theta_src)  -  np.cos(theta_stn)*np.sin(theta_src)*np.cos(phi_stn - phi_src))

    rot2 = (1/np.sin(dist)) * (np.sin(theta_src)*np.sin(phi_stn - phi_src))

    # Conversion from RTP --> ZNE (where R=Z, 2D rotation) appears to use the following matrix:
    #   [N, E]' = [-rot1, -rot2; -rot2, rot1][T, P]' where T and P are theta, Phi
    #   Below we shall name the rotation matrix Q:
    # Hence to get the T and P matrix we should be multiplying [N,E] by the inverse of Q:
    Q    = np.array([[-rot1, -rot2], [rot2, -rot1]])
    Qinv = np.linalg.inv(Q)


    if method == "NE->TP":
        data_NE = np.array([traceN, traceE])
        data_TP = np.matmul(Qinv, data_NE)

        # Now writing back to stream:
        old_chls = ["N", "E"]
        new_chls = ["T", "P"]
        for i in range(2):
            if new_chls[i] == "P":
                data_TP[i, :] = data_TP[i,:]*(-1)

    else:
        raise ValueError("Currently method must be NE->TP")

    return data_TP[0,:], data_TP[1,:]