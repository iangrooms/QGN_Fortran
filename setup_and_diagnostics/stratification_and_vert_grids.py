import numpy as np
import scipy.integrate as integrate
import scipy.optimize as optimize
import scipy.linalg as linalg

def cauchy_dist(z,z0,w):
    """ cauchy distribution helper for defining synthetic ocean stratification """
    return (w/np.pi)/((z-z0)**2+w**2)

def get_rho(Htot,N2,zs,f=1e-4,rho_top=1025.5,rho0=1030):
    """
    Integrate to obtain densitites for stratification profile densities;
    
    Parameters:
    N2: buoyancy frequency function
    zs: geopotential heights of cell interfaces at which to evaluate density
    rho_top: density at the surface
    rho0: reference density

    Return:
    Densities at provided heights zs
    """
    
    import scipy.integrate as integrate
    
    rho = np.zeros(len(zs))
    for i in range(0,len(zs)):
        rho[i],_ =  integrate.quad(N2,zs[i],1)
        
    g = 9.81
    return (f**2*rho0*Htot/g)*rho + rho_top

def get_alpha_grid(nz,alpha,N2):
    """
    return an equispaced vertical grid with nz layers in the given alpha coordinate described in :
    Rachel Robey and Ian Grooms. Continuous and discrete baroclinic modes in continuously varying stratification.
    SIAM Journal on Applied Mathematics, 84(6):2502–2521, 2024. https://doi.org/10.1137/24M1648181.

    Parameters:
    nz: number of layers
    alpha: vertical coordinate parameter, any float in [0,2]. In particular:
           alpha = 0 uses geopotential
           alpha = 1 uses Charney coordinate
           alpha = 2 uses isopycnal coordinate

    Return
    z: cell interfaces of the grid
    zc: cell centers
    h: cell heights
    """
    
    SS = lambda t: np.sqrt(N2(t)) ** alpha
    normalization_constant = integrate.quad(SS,0,1)[0]
    z = np.zeros(nz+1)
    z[nz] = 1
    for i in range(1,nz):
        z[i] = optimize.fsolve(lambda t: integrate.quad(SS,0,t)[0] - normalization_constant*i/nz,i/nz)

    return z

# Half Chebyshev grid in the Charney coordinate
def get_CC_grid(nz,alpha,N2):
    """
    
    return Chebychev Tro, S., Robey, R., and Grooms, I. Ekman-driven buoyancy flux in quasigeostrophic flow. Submitted to J. Fluid Mech., Apr 2025.
    """
    
    SS = lambda t: np.sqrt(N2(t)) ** alpha
    normalization_constant = integrate.quad(SS,0,1)[0]
    
    sc = 1 - np.cos(np.pi * np.linspace(0,1,nz+1) / 2) # Interfaces at the interior points of a half-Chebyshev grid
    
    z = np.zeros(nz+1)
    z[nz] = 1
    for i in range(1,nz):
        z[i] = optimize.fsolve(lambda t: integrate.quad(SS,0,t)[0] - normalization_constant*sc[i],i/nz)
        
    return z


def get_modes(z,h,N2strat):
    """
    return baroclinic radii (eigenvalues), baroclinic modes (eigenvectors), and stretching matrix L

    Parameters:
    N2strat: function of the stratification/buoyancy frequency with respect to depth
    z: depth of grid interfaces
    h: heights of cells

    Return:
    vals: eigenvalues of the stratified stretching matrix
    vecs: eigenvectors of the stretching matrix (i.e. baroclinic modes)
    L: the stretching matrix
    """
    
    nz = len(h)
    assert len(z) == nz+1
    
    DS = np.diag(1/N2strat(z[1:-1]))
    D0 = np.zeros((nz-1,nz))
    D0[:,:-1] = 2*np.diag(1/(h[:-1]+h[1:]))
    D0[:,1:] = D0[:,1:] - 2*np.diag(1/(h[:-1]+h[1:]))
    D = DS@D0
    L = (1/h)*(np.vstack((np.zeros((1,nz)),D))-np.vstack((D,np.zeros((1,nz)))))
    L = L.T # Use C/Python's indexing convention rather than Fortran/Matlab.
    
    (vals,vecs) = np.linalg.eig(L)
    
    ## sort and normalize modes / eigenvalues
    ind = np.argsort(-vals)
    vecs = vecs[:,ind]
    vals = vals[ind]
    
    for k in range(0,nz):
        tmp = np.sum(h*vecs[:,k]**2)
        vecs[:,k] = vecs[:,k] / np.sqrt(tmp)
        
    for k in range(0,nz):
        vecs[:,k] = np.sign(vecs[0,k])*vecs[:,k]
        
    return vals, vecs, L

def get_zeros(f,x):

    i0 = np.array(np.where(f == 0))[0]
    ix = np.array(np.where(f[1:]*f[:-1] < 0))[0]
    diff = np.abs(f[1:] - f[:-1])

    w = np.abs(f[ix])/diff[ix]

    i0 = np.sort(np.concatenate([i0,x[ix]*(1-w)+x[ix+1]*w]))
    return i0

def get_roots_grid(nz,N2):
    """
    return vertical grid with interfaces placed at the roots of the
    oscillations of the second largest (nz-1) baroclinic mode

    Parameters:
    nz: number of layers

    Return:
    z: cell interfaces of the grid
    zc: cell centers
    h: cell heights
    """

    zi_ref = np.linspace(0,1,1025)
    zc_ref = (zi_ref[1:]+zi_ref[:-1])/2.
    h_ref = np.ones(1024)*(zi_ref[1]-zi_ref[0])

    _, vecs, _ = get_modes(zi_ref,h_ref,N2)

    nth_mode = vecs[:,nz-1]
    zi = np.flip(np.concatenate([[0],get_zeros(nth_mode,zc_ref),[1]]))
    
    return zi

def get_extrema_grid(nz,N2):
    """
    return vertical grid with interfaces placed at the extrema
    of the oscillations of the largest (nz) baroclinic mode

    Parameters:
    nz: number of layers

    Return:
    z: cell interfaces of the grid
    zc: cell centers
    h: cell heights
    """
    
    from scipy.signal import argrelextrema

    zi_ref = np.linspace(0,1,1025)
    zc_ref = (zi_ref[1:]+zi_ref[:-1])/2.
    h_ref = np.ones(1024)*(zi_ref[1]-zi_ref[0])

    _, vecs, _ = get_modes(zi_ref,h_ref,N2)

    nth_mode = vecs[:,nz]
    iex = argrelextrema(np.abs(nth_mode), np.greater)
    zi = np.flip(np.concatenate([[0],zi_ref[iex],[1]]))
    
    return zi

def get_stewart_grid(Htot,H=6e3, dzd = 199.1, min_dz=2.3, depfac=1.01):

    """
    Grid generation from 
    Stewart, Kial & Hogg, A.McC & Griffies, Stephen & Heerdegen, A.P. & Ward, M.L. & Spence, P. & England, Matthew. (2017).
    Vertical resolution of baroclinic modes in global ocean models. Ocean Modelling. 113. 10.1016/j.ocemod.2017.03.012.
    as posted at
    https://github.com/kialstewart/vertical_grid_for_ocean_models

    Parameters:
    H: maximum depth of your ocean (approximately)? in meters
    dzd: maximum grid spacing (the grid spacing at the deepest point in the ocean) in meters
    min_dz: minimum grid spacing (the grid spacing at the ocean surface) in meters
    depfac: tune sharpness of the hyperbolic tangent (<1 is sharp, 1 is neutral, >1 is gentle) / total number of levels

    """

    ################
    # start the build
    ################
    
    # import netCDF4 as nc
    # import numpy as np
    
    # define the functional form of the vertical grid
    epsilon = 0.001 # this is a small number needed to begin the iteration
    def f_all(kk):
        return np.tanh(np.pi*((kk)/(H*depfac)))*(dzd)+epsilon # the function is {tanh(pi*m/H)*dz_max + epsilon}, which is epsilon at the surface and dz_max at H
    
    # make the first two entries of the initial grid; these will be 0 and epsilon for both z and dz
    delta_z = [0,epsilon*1.0]
    prop_z = [0,epsilon*1.0]
    
    # this is where the magic happens: an iterative process that takes a step from the current end (deepest point) of the grid along the function to find the next point
    while prop_z[-1]+delta_z[-1] < 1.2*H:
        aa = np.linspace(1.0,1.5,10000)
        bb = np.zeros([len(aa)])
        loopkill = 1.0
        ii = 0
        while loopkill > 0:
            bb[ii] = (f_all(prop_z[-1]+(delta_z[-1]*aa[ii])))-(delta_z[-1]*aa[ii])
            loopkill = bb[ii]
            ii += 1
        aa_bb = np.polyfit(aa[:ii-1],bb[:ii-1],1)
        dznew = (delta_z[-1]*(np.abs(aa_bb[1]/aa_bb[0])))
        delta_z = np.append(delta_z,dznew)
        prop_z = np.append(prop_z,(prop_z[-1]+delta_z[-1]))
    
    # now that we have an initial grid that follows the desired functional form we need to relocate it vertically so that the grid spacing at the surface is min_dz
    new_surf = np.max(np.where(delta_z<min_dz)) # find where the initial grid is min_dz (the surface resolution)
    real_prop_z = prop_z[new_surf:]-prop_z[new_surf] # make a new grid that shifts the initial grid vertically
    # real_delta_z = delta_z[new_surf:] # make a new dz for this new grid
    real_prop_z = real_prop_z[np.where(real_prop_z<H)] # cut the new grid off at desired depth, H
    # real_delta_z = real_delta_z[np.where(real_prop_z<H)] # and the new dz too
    
    # print("SUCCESS!! Created vertical grid of", len(real_prop_z), "levels with grid spacing from %.2f m to %.2f m"%(real_delta_z[0],real_delta_z[-1]))
    
    zi_stwt = (1-np.concatenate([real_prop_z[np.where(real_prop_z < Htot)]/Htot,[1]]))
    
    return zi_stwt

def get_MOM6_zGrid(Htot):

    """
    A geopotential grid from MOM6; used in 
    Marques, G. M., Shao, A. E., Bachman, S. D., Danabasoglu, G., & Bryan, F. O. (2023). Representing eddy diffusion in the surface
    boundary layer of ocean models with general vertical coordinates. Journal of Advances in Modeling Earth Systems, 15, e2023MS003751.
    https://doi.org/10.1029/2023MS003751 
    """

    h_mom6 = np.array([2.5, 2.5, 2.5, 2.5, 2.77, 3.38, 4.01, 4.65, 5.29, 5.95, 6.61, 7.28, 7.97, 8.66, 9.37,\
                      10.08, 10.81, 11.54, 12.29, 13.06, 13.85, 14.69, 15.59, 16.56, 17.61, 18.76, 20.02,\
                      21.42, 23.0, 24.77, 26.79, 29.1, 31.76, 34.87, 38.5, 42.79, 47.9, 54.01, 61.37, 70.25,\
                      80.95, 93.75, 108.8, 126.04, 145.04, 164.81, 184.05, 201.34, 215.66, 226.64, 234.5,\
                      239.84, 243.31, 245.52, 246.88, 247.72, 248.23, 248.54, 248.73, 248.84, 248.64, 248.68,\
                      248.71, 248.72, 248.73])

    h_mom6 = (h_mom6[np.where(np.cumsum(h_mom6)<Htot)])/Htot
    zi_mom6 = (1-np.concatenate([[0],np.cumsum(h_mom6),[1]]))

    return zi_mom6

def get_OM4_isoGrid(Htot, N2, f=1e-4, beta=2e-11, rho_top=1025., shft=9):
    """
    Adcroft, A., Anderson, W., Balaji, V., Blanton, C., Bushuk, M., Dufour, C. O., et al. (2019). The GFDL global ocean and
    sea ice model OM4.0: Model description and simulation features. Journal of Advances in Modeling Earth Systems, 11, 3167–3211.
    https://doi.org/10.1029/2019MS001726 
    
    Grid using isopycnal portion of the OM4 hybrid coordinate. A tuneable shift was introduced to improve the coverage of the target densities over 
    the range of values in the suppliedcover the density range for the supplied synthetic ocean stratification. Returns geopotential coordinates 
    determined by the intersections of the density profile with the target values given.
    
    The target potential density values are referenced to 2,000 dbar, $\rho_2$ with 
    $$\rho_p=\rho(S,\theta,p_r+0.01(p-p_r))$$

    Parameters:
    f: local Coriolis
    beta: local Coriolis gradient
    rho_top: density at the surface
    tune_shift: tunable parameter to tweak target range
    """
    
    rhoiso = np.array([1010.0, 1014.3034, 1017.8088, 1020.843, 1023.5566, 1025.813, 1027.0275, 1027.9114, 1028.6422,\
                   1029.2795, 1029.852, 1030.3762, 1030.8626, 1031.3183, 1031.7486, 1032.1572, 1032.5471,\
                   1032.9207, 1033.2798, 1033.6261, 1033.9608, 1034.2519, 1034.4817, 1034.6774, 1034.8508,\
                   1035.0082, 1035.1533, 1035.2886, 1035.4159, 1035.5364, 1035.6511, 1035.7608, 1035.8661,\
                   1035.9675, 1036.0645, 1036.1554, 1036.2411, 1036.3223, 1036.3998, 1036.4739, 1036.5451,\
                   1036.6137, 1036.68, 1036.7441, 1036.8062, 1036.8526, 1036.8874, 1036.9164, 1036.9418,\
                   1036.964, 1036.9857, 1037.0052, 1037.0236, 1037.0409, 1037.0574, 1037.0738, 1037.0902,\
                   1037.1066, 1037.123, 1037.1394, 1037.1558, 1037.1722, 1037.1887, 1037.206, 1037.2241,\
                   1037.2435, 1037.2642, 1037.2866, 1037.3112, 1037.3389, 1037.3713, 1037.4118, 1037.475,\
                   1037.6332, 1037.8104, 1038.0])

    zs = np.linspace(0,1,1024)
    rho = get_rho(Htot,N2,zs,rho_top=rho_top)
    
    ni = np.where( (rhoiso-shft >= rho[-1]) & (rhoiso-shft <= rho[0]) )[0]
    intersections = np.zeros(len(ni))
    for i in range(0,len(ni)):
        val = (rhoiso[ni[i]] - shft - rho_top)*9.81/(f**2*1030*Htot)
        intersections[i] = optimize.fsolve(lambda t: integrate.quad(N2,t,1)[0] - val,1-i/len(ni))
    # ni = intersections[np.where(intersections > 0)]
    
    zi_iso = np.concatenate([[1],intersections,[0]])

    return zi_iso

