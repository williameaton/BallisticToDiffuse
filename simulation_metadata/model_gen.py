# File script to generate scattering models for synthetic studies. Could be converted into one main function but easier
# to run as a file script from the command line.
# All classes and functions are defined at top
# Parameters can be defined at the bottom of script for model. 

# Imports
import numpy as np
import math
import netCDF4 as nc

# ----------------------------------------------------------------------------------------------------------------------

# Define classes


# Class for creating 3D model objects:
class model(object):
    def __init__(self, x_lim, y_lim, z_lim, elements_per_wavelength, dominant_freq, min_velocity, oversaturation=1):
        # Note that oversaturation can be used to define a mesh resolution above that which is necessary from the
        # freq/vel etc.
        self.x_lim = x_lim  # x-axis limit of model 
        self.y_lim = y_lim  # y-axis limit of model  
        self.z_lim = z_lim  # z-axis limit of model 
        self.elements_per_wavelength = elements_per_wavelength # epw
        self.freq = dominant_freq  # dominant frequency for the model which defines the mesh spacing 
        self.min_velocity = min_velocity # minimum velocity also defines the grid spacing 

        # Calculate min wavelength from frequency and min velocity
        self.min_wavelength = min_velocity/dominant_freq

        # Calculate the number of elements in each dimension based on domain size, min. wavelength and elements per
        # wavelength
        self.nx = oversaturation*math.ceil((x_lim[1] - x_lim[0])*elements_per_wavelength/self.min_wavelength)
        self.ny = oversaturation*math.ceil((y_lim[1] - y_lim[0])*elements_per_wavelength/self.min_wavelength)
        self.nz = oversaturation*math.ceil((z_lim[1] - z_lim[0])*elements_per_wavelength/self.min_wavelength)

        # Store data of initial model parameters
        self.padding = np.array([0, 0, 0])
        self.original_shape = np.array([self.nx, self.ny, self.nz])


        # Calculate the length of each dimension:
        self.x_length = self.x_lim[1] - self.x_lim[0]
        self.y_length = self.y_lim[1] - self.y_lim[0]
        self.z_length = self.z_lim[1] - self.z_lim[0]

        # Calculate spatial step sizes
        self.dx = self.x_length / self.nx
        self.dy = self.y_length / self.ny
        self.dz = self.z_length / self.nz

        # Define a default background 3D model arrays for Rho, Vp and Vs which are zeros 3D arrays of the correct size:
        self.bm_rho = np.zeros((self.nx, self.ny, self.nz))
        self.bm_vp  = np.zeros((self.nx, self.ny, self.nz))
        self.bm_vs  = np.zeros((self.nx, self.ny, self.nz))

        self.unpadded_n = np.array([self.nx, self.ny, self.nz])


# Class for generating sphere objects which may then be injected into the model 
class sphere():
    def __init__(self, sphere, vp, vs, rho, radius):
        self.sphere = sphere # Array that holds the sphere (3D array of 1 or 0 entries that has spherical structure)
        self.vp = vp # Vp for sphere
        self.vs = vs # Vs for sphere
        self.rho = rho # Density for sphere
        self.centre = np.array([0,0,0]) # coordinates for the centre of sphere within the model
        self.n_centre = np.array([0,0,0]) 
        self.radius = radius # Sphere radius 

        # Get the index values of the sphere array for the centre point of that array
        self.sa_centre = np.array([0,0,0])
        for i in range(3):
            self.sa_centre[i] = np.floor(np.asarray(self.sphere.shape)[i]/2)


    # Functions to update sphere properties: 
    def set_radius(self, radius):
        self.radius = radius

    def update_vp(self, new_vp):
        self.vp = new_vp

    def update_vs(self, new_vs):
        self.vs = new_vs

    def update_rho(self, new_rho):
        self.rho = new_rho

    def set_centre(self, centre, model, print_conf='n'):
        self.centre = centre
        
        # Initialise centre:
        self.n_centre = np.array([0, 0, 0])

        self.n_centre[0] = model.unpadded_n[0] * centre[0] // model.x_length
        self.n_centre[1] = model.unpadded_n[1] * centre[1] // model.y_length
        self.n_centre[2] = model.unpadded_n[2] * centre[2] // model.z_length

        if print_conf.upper() == "Y":
            print("Centre of sphere set at", centre, "with normalised indices", self.n_centre)

    def _update_sph_centre_index(self, new_index):
        self.sa_centre = new_index

# ----------------------------------------------------------------------------------------------------------------------
# ______________________________________________________________________________________________________________________
# Define functions:
# ______________________________________________________________________________________________________________________

def gen_sphere(m, radius, RHO, VP, VS, print_time='n'):
    # INPUT VARIABLES: 
    # m [Model object]      - Instance of Model class required so grid spacing data etc can be extracted
    # radius [int, float]   - Radius of sphere
    # RHO, VP, VS [float]   - Physical properties of sphere 
    
    # OUTPUT: 
    # sphere_instance       - Returns instance of sphere object
    
    
    # Calculate domain sizes for scaling:
    x_len = m.x_lim[1] - m.x_lim[0]
    y_len = m.y_lim[1] - m.y_lim[0]
    z_len = m.z_lim[1] - m.z_lim[0]

    # Calculate spatial step sizes
    dx = x_len / m.nx
    dy = y_len / m.ny
    dz = z_len / m.nz

    r_original = radius

    # Calculate the number of iterations based on radius and element size:
    x_loop = int(radius // dx)
    y_loop = int(radius // dy)
    z_loop = int(radius // dz)

    # Create spare array that holds sphere info for duplicate spheres
    # This will store values of '1' or '0' for whether the element is within the sphere radius (1) or not (0)
    sph = np.zeros((int(x_loop+1), int(y_loop+1), int(z_loop+1)))


    # Calculating the valid array elements for a circle in the positive x,y,z octet
    # Loop over all dimensions:
    for k in np.arange(0, z_loop + 1):
        for j in np.arange(0, y_loop + 1):
            for i in np.arange(0, x_loop + 1):

                # check element is within radius:
                element_radius = (((i * dx) ** 2) + ((j * dy) ** 2) + ((k * dz) ** 2)) ** 0.5

                if element_radius <= radius:

                    sph[i, j, k] = 1

    # This loop gives us the positive octet - now need to mirror for all other octets:
    # This reflects s along the first column and row as mirror axes and adds sph_x - 1 etc elements (ie ignoring the
    # row/column being used as the mirror axis)
    # Initialise matrix
    sph_full = np.zeros((2*int(x_loop+1)-1, 2*int(y_loop+1)-1, 2*int(z_loop+1)-1))
    sph_full = np.lib.pad(sph, ((x_loop, 0), (y_loop, 0), (z_loop, 0)), 'reflect')

    # Create instance of the sphere class
    sphere_instance = sphere(sph_full, VP, VS, RHO, radius)
    return sphere_instance


# Alternative sphere generation function in which an extra spherical shell is added to the outer of the original sphere
# with a gaussian taper (i.e. perturbation strength tapers from 100 % to 0 % of its original value) within this extra
# sphere.
def gen_sphere_gaussian(m, radius, RHO, VP, VS, buffer=3, falloff=1, print_time='n'):
    # INPUT VARIABLES: 
    # m [Model object]      - Instance of Model class required so grid spacing data etc can be extracted
    # radius [int, float]   - Radius of sphere
    # RHO, VP, VS [float]   - Physical properties of sphere 
    # buffer [int]          - Extra shell size (buffer) 
    # falloff [int]         - Falloff of gaussian taper (i.e related to the std dev of the gaussian)
    # print_time [str]      - Used for develop to test efficiency of functions - use 'y' or 'yes' to print
    
    # OUTPUT: 
    # sphere_instance       - Returns instance of sphere object
    
    # Calculate domain sizes for scaling:
    x_len = m.x_lim[1] - m.x_lim[0]
    y_len = m.y_lim[1] - m.y_lim[0]
    z_len = m.z_lim[1] - m.z_lim[0]

    # Calculate spatial step sizes
    dx = x_len / m.nx
    dy = y_len / m.ny
    dz = z_len / m.nz

    r_original = radius
    radius += falloff*dx*buffer

    # Calculate the number of iterations based on radius and element size:
    x_loop = int(radius // dx)
    y_loop = int(radius // dy)
    z_loop = int(radius // dz)

    # Create spare array that holds sphere info for duplicate spheres
    # This will store values of '1' or '0' for whether the element is within the sphere radius
    sph = np.zeros((int(x_loop+1), int(y_loop+1), int(z_loop+1)))

    # calculating values for gaussian:
    gauss_sd = dx

    # Calculating the valid array elements for a circle in the positive x,y,z octet
    # Loop over all dimensions:
    for k in np.arange(0, z_loop + 1):
        for j in np.arange(0, y_loop + 1):
            for i in np.arange(0, x_loop + 1):

                # check element is within radius:
                element_radius = (((i * dx) ** 2) + ((j * dy) ** 2) + ((k * dz) ** 2)) ** 0.5

                if element_radius <= radius:
                    # Apply parameter to the element
                    if element_radius <= r_original:
                        sph[i, j, k] = 1
                    else:
                        sph[i,j,k] = 1 * np.exp(-0.5 * (element_radius - r_original) ** 2 / gauss_sd ** 2)
    # This loop gives us the positive octet - now need to mirror for all other octets:
    # This reflects s along the first column and row as mirror axes and adds sph_x - 1 etc elements (ie ignoring the
    # row/column being used as the mirror axis)
    # Initialise matrix
    sph_full = np.zeros((2*int(x_loop+1)-1, 2*int(y_loop+1)-1, 2*int(z_loop+1)-1))
    sph_full = np.lib.pad(sph, ((x_loop, 0), (y_loop, 0), (z_loop, 0)), 'reflect')

    # Create instance of the sphere class
    sphere_instance = sphere(sph_full, VP, VS, RHO, radius)
    return sphere_instance


    # Get non-zero values of injects:
    truth_vp = vp_inj != 0
    truth_vs = vs_inj != 0
    truth_rho = rho_inj != 0

    model.bm_vp[lb[0]:ub[0], lb[1]:ub[1] ,lb[2]:ub[2]][truth_vp] = vp_inj[truth_vp]
    model.bm_vs[lb[0]:ub[0], lb[1]:ub[1] ,lb[2]:ub[2]][truth_vs ]= vs_inj[truth_vs]
    model.bm_rho[lb[0]:ub[0], lb[1]:ub[1] ,lb[2]:ub[2]][truth_rho] += rho_inj[truth_vp]


# ______________________________________________________________________________________________________________________
def centre_create(model, mfp, domain, spec_domain,  type='edge'):

    if domain.upper() == 'SINGLE':
        X, Y, Z = np.mgrid[model.x_lim[0]:model.x_lim[1]+0.1:mfp,
                           model.y_lim[0]:model.y_lim[1]+0.1:mfp,
                           model.z_lim[0]:model.z_lim[1]+0.1:mfp]

        # Centre in x and y:
        x_centre_max = X[-1,0,0]
        y_centre_max = Y[0,-1,0]
        z_centre_max = Z[0,0,-1]

        x_add = (model.x_lim[1]-x_centre_max)/2
        y_add = (model.y_lim[1]-y_centre_max)/2
        z_add = (model.z_lim[1]-z_centre_max)/2

        X += x_add
        Y += y_add
        Z += z_add

        xyz = np.vstack((X.flatten(), Y.flatten(), Z.flatten())).T
        return(xyz)
# _____________________________________________________________________________________________________________________________________________________


def slice_sphere(sph, model, print_conf='N'):
    # Extract the sphere array
    sph_array = sph.sphere
    # extract sphere centre:
    n_centre = sph.n_centre
    # Calculate the dimensions of sphere array:
    sph_dim = np.asarray(sph_array.shape)
    # Get the index for the centre of sphere array - may need updating if slice elemets with lower indices than it
    sa_centre = sph.sa_centre
    # Get dimensions of model array:
    mod_dim = [model.nx, model.ny, model.nz]
    # Calculating the number of elements either side of the centre in the sphere array in each dimension:
    el_count = np.floor(sph_dim/2)

    # For the lower bounds:
    # Calculating difference between the number of elements left/beneath the centre in the sph_array and the index of
    # the normalised centre position
    lb_out = np.floor(el_count - n_centre)
    # in x
    if lb_out[0] > 0:
        sph_array = sph_array[int(lb_out[0]):, :, :]
        sa_centre[0] = sa_centre[0] - lb_out[0]
    # in y
    if lb_out[1] > 0:
        sph_array = sph_array[:, int(lb_out[1]):, :]
        sa_centre[1] = sa_centre[1] - lb_out[1]

    # in z
    if lb_out[2] > 0:
        sph_array = sph_array[:, :, int(lb_out[2]):]
        sa_centre[2] = sa_centre[2]-lb_out[2]

    # For the upper bound:
    # If this is below zero then it means there are more on one side of the sph_array than there are elements between
    # the centre point and the edge of the domain - need slicing from RHS
    ub_out = np.floor(mod_dim - (n_centre + 1) - el_count)
    # in x
    if ub_out[0] < 0:
        sph_array = sph_array[:int(ub_out[0]), :, :]
    # in y
    if ub_out[1] < 0:
        sph_array = sph_array[:, :int(ub_out[1]), :]
    # in z
    if ub_out[2] < 0:
        sph_array = sph_array[:, :, :int(ub_out[2])]

    # Create sliced sphere instance:
    sphere_sliced = sphere(sph_array, vp=sph.vp, vs=sph.vs, rho=sph.rho, radius=sph.radius)
    # Need to set centre to be the same as for the normal sphere class:
    sphere_sliced.set_centre(centre=sph.centre, model=model )

    # Update the index location of the sphere centre within the sphere array:
    sphere_sliced._update_sph_centre_index(sa_centre)

    return sphere_sliced
# ______________________________________________________________________________________________________________________
def inject_sphere(model, sphere, print_conf="N"):
    # Lower bound index is given by the index location of the sphere centre (in domain coordinates) - index location of
    # the sphere's centre within the sliced sphere array
    lb = sphere.n_centre - sphere.sa_centre # array for [z,y,x]
    # Upper bound is given by centre of the sphere (in domain) + shape of the sliced array - index loc of the sphere
    ub = sphere.n_centre + np.asarray(sphere.sphere.shape) - sphere.sa_centre
    # Creating array to inject for each parameter:
    vp_inj = sphere.sphere * sphere.vp
    vs_inj = sphere.sphere * sphere.vs
    rho_inj = sphere.sphere * sphere.rho

    non_zero = vp_inj != 0
    # inject:
    model.bm_vp[lb[0]:ub[0], lb[1]:ub[1] ,lb[2]:ub[2]] = vp_inj
    model.bm_vs[lb[0]:ub[0], lb[1]:ub[1] ,lb[2]:ub[2]] = vs_inj
    model.bm_rho[lb[0]:ub[0], lb[1]:ub[1] ,lb[2]:ub[2]] = rho_inj

    return model
# ______________________________________________________________________________________________________________________
def addSphere(sph, model, print_time='n'):

    # Extract sph_input:
    sph_input = sphere(sph.sphere, sph.vp, sph.vs, sph.rho, sph.radius)
    # update sph_input centre:
    sph_input.set_centre(sph.centre, model)

    # Slice the inputted sphere to the correct array size
    sphere_sliced = slice_sphere(sph_input, model)

    # Inject into model:
    final = inject_sphere(model, sphere_sliced)


    return final
# ______________________________________________________________________________________________________________________
def spaced_spheres(sphere, model, mfl, print_time='n', domains='SINGLE', spec_domain=[], type='edge'):

    # Create an array holding the location of sphere centres:
    sc = centre_create(model, mfl, domains, spec_domain=spec_domain,  type=type, )

    for i in range(sc.shape[0]):
        # print("Adding sphere at ", sc[i, :])
        # Set centre of the sphere:
        sphere.set_centre(sc[i, :], model, print_conf='n')
        # Add sphere to model:
        final = addSphere(sphere, model, print_time='n')

    print("Added", sc.shape[0], "spheres to model")
    return final
# ______________________________________________________________________________________________________________________
def writeNetCDF(m, filename):

    f = nc.Dataset(filename, 'w', format='NETCDF4')
    # Create dimension arrays
    x_array = np.linspace(-m.x_lim[1]/2, m.x_lim[1]/2, m.nx)
    y_array = np.linspace(-m.y_lim[1]/2, m.y_lim[1]/2, m.ny)
    z_array = np.linspace(-m.z_lim[0], m.z_lim[1], m.nz)
    # We now create the dimensions:
    x_dim = f.createDimension('x_dim', m.nx)
    y_dim = f.createDimension('y_dim', m.ny)
    z_dim = f.createDimension('z_dim', m.nz)
    # Creating the variables:
    x   = f.createVariable('x', 'f4', ('x_dim',))
    y   = f.createVariable('y', 'f4', ('y_dim',))
    z   = f.createVariable('z', 'f4', ('z_dim',))
    v_rho = f.createVariable('rho', 'f4', ('x_dim', 'y_dim', 'z_dim',))
    v_vp  = f.createVariable('vp', 'f4', ('x_dim', 'y_dim', 'z_dim',))
    v_vs  = f.createVariable('vs', 'f4', ('x_dim', 'y_dim', 'z_dim',))
    # Assigning values to the variables:
    x[:] = x_array
    y[:] = y_array
    z[:] = z_array
    v_rho[:,:,:] = m.bm_rho
    v_vp[:,:,:]  =  m.bm_vp
    v_vs[:,:,:]  =  m.bm_vs
    print('Data written to file ', filename)
    f.close()
# ______________________________________________________________________________________________________________________

# EXAMPLE INPUT PARAMETERS TO GENERATE MODEL 
# Would reccommend using SI units only - e.g below uses 250 km for x and y with 125 km for z dimension

# RUN SCRIPT:
# Model dimensions
x = np.array([0, 100000])
y = np.array([0, 100000])
z = np.array([0, 100000])
# Model parameters

freq = 2       # Hz
min_vel = 3000 # m/s
epw = 3        # Note that changing the EPW will dramatically change the size of the output model file


# Create model class
m = model(x,y, z, epw, freq, min_vel, oversaturation=1)

# Define perturbation
ptb = -0.2
vp_ptb = -0.2
rho_ptb = -0.2


sphere_rad_wavelengths = 1 # Sphere radius in wavelengths
mfp = 2                    # In wavelengths (mfp is distance between spheres - mean free path)


# Generate spheres:
sph = gen_sphere(m, sphere_rad_wavelengths*m.min_wavelength, rho_ptb, vp_ptb, ptb)
out = spaced_spheres(sph, m, mfp*m.min_wavelength, print_time='n')

scr = f'p{str(ptb)}_{str(freq)}hz_{str(mfp)}_mfp_{str(sphere_rad_wavelengths)}_rad.nc'

writeNetCDF(out, scr)
