import re
import os
import glob
import struct

import numpy as np
import xarray as xr
import scipy.linalg as linalg
from scipy.fft import fft2, fftfreq

def get_timesteps(case_dir,logfile='out.txt'):
    """
    Get the size of the timesteps from the raw outfile from QGN

    Input
    case_dir :: path to simulation directory including raw QGN log file, default out.txt
    logfile :: name of the output log file from the simulation, default 'out.txt'

    Output
    time :: total simulation time for each logged iteration in days (following QGN output convention)
    tstep :: timestep for each logged iteration in seconds (following QGN output convention)
    """
    
    import re
    
    match_number = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?')

    time = [0]       ## total simulation time from start in days
    tstep = [np.nan] ## time step in seconds
    with open(case_dir+logfile) as f:
        l = 0
        while True:
            line = f.readline()
            if not line:
                return tstep
            l = l+1
            if 'Time since inception' in line:
                break
                
        while True:
            res = [float(x) for x in re.findall(match_number, line)]
            try:
                time.append(res[0])
                tstep.append(res[1])
            except:
                break
            
            while True:
                line = f.readline()
                if not line:
                    return time, tstep
                l = l+1
                if 'Time since inception' in line:
                    break
                
    return time, tstep

def get_mode_order(case_dir,nz,logfile='out.txt'):
    """
    Scrape the eigenvalues/modes from the QGN output file to get the order to sort the modal output

    Input
    case_dir :: path to simulation directory including raw QGN log file, default out.txt
    nz :: number of vertical layers
    logfile :: name of the output log file from the simulation, default 'out.txt'

    Output
    order :: an indexing order of length nz
    """
    
    ## scrape outfile for eigenvalues / mode ordering
    eigvals = []
    eigflag = False
    with open(case_dir+logfile) as f:
        while len(eigvals) < nz:
            line = f.readline()
            if eigflag:
                if "--------" in line:
                    break
                else:
                    eigvals.extend(line.split())

            if "The eigenvalues are" in line:
                eigflag = True

    eigvals = np.array(eigvals,dtype=float)
    order = np.argsort(eigvals)[::-1]
    
    assert len(order) == nz

    return order

def get_output(case_dir,it,Htot,x,y):
    """
    Read QGN output into an xarray dataset. Directory needs the following:
    zc.dat
    H.dat
    q.%09d.dat (iteration stamp formatted in QGN)
    p.%09d.dat (iteration stamp formatted in QGN)


    Input
    case_dir :: path to run directory including listed input and output files for the simulation
    it :: iteration number
    Htot :: height of the simulation domain for scaling grid point depths
    x, y :: array of horizontal coordinates in x and y directions respectively
    nz :: 


    Output
    ds :: xarray of q and psi output with vertical and horizontal coordinates
    """
    
    nx = len(x)
    ny = len(y)

    H = np.fromfile(case_dir+'H.dat')
    
    zc = (np.fromfile(case_dir+'zc.dat')-1)*Htot
    zi = (-np.cumsum(np.concatenate([[0],H])))*Htot

    nz = len(zc)

    q = np.reshape(np.fromfile(case_dir+'q.%09d.dat'%it),(nz,ny,nx))
    p = np.reshape(np.fromfile(case_dir+'p.%09d.dat'%it),(nz,ny,nx))
        
    ds = xr.Dataset({'q':(['zc','y','x'],q),\
                     'p':(['zc','y','x'],p),},\
                     coords={'zc':zc,'zi':zi,'x':x,'y':x})

    return ds

def get_modal_output(case_dir,it,x,y,nz,logfile='out.txt'):
    """
    Read QGN modal output into an xarray dataset. Directory needs the following:
    q_mode.%09d.dat (iteration stamp formatted in QGN)
    p_mode.%09d.dat (iteration stamp formatted in QGN)
    J_mode.%09d.dat (iteration stamp formatted in QGN)


    Input
    case_dir :: path to run directory including listed input and output files for the simulation
    it :: iteration number
    x, y :: array of horizontal coordinates in x and y directions respectively


    Output
    ds :: xarray of q and psi and J[psi,q] output in vertical modes with vertical and horizontal coordinates
    """
    nx = len(x)
    ny = len(y)
    
    order = get_mode_order(case_dir,nz,logfile=logfile)
    
    qmode = np.reshape(np.fromfile(case_dir+'q_mode.%09d.dat'%it),(nz,ny,nx))[order,:,:]
    pmode = np.reshape(np.fromfile(case_dir+'p_mode.%09d.dat'%it),(nz,ny,nx))[order,:,:]
    Jmode = np.reshape(np.fromfile(case_dir+'J_mode.%09d.dat'%it),(nz,ny,nx))[order,:,:]
        
    ds = xr.Dataset({'qmode':(['m','y','x'],qmode),\
                     'pmode':(['m','y','x'],pmode),\
                     'Jmode':(['m','y','x'],Jmode)},\
                     coords={'m':np.arange(nz),'x':x,'y':x})

    return ds

def add_prof_vars(case_dir,ds,Htot):
    """
    Add stratification and layer height to xarray dataset. Directory needs the following:
    S.dat
    H.dat


    Input
    case_dir :: path to run directory including listed input files for the simulation
    ds :: the dataset to add the variables to
    Htot :: height of the simulation domain for scaling grid layer heights

    Output
    ds :: xarray with added S and h variables
    """

    S = np.fromfile(case_dir+'S.dat')
    H = np.fromfile(case_dir+'H.dat')*Htot
    
    ds['S'] = xr.DataArray(S,dims=['zi'])
    ds['h'] = xr.DataArray(H,dims=['zc'])

    return ds

def get_KE(p,Lx=1,Ly=1):
    """
    At a fixed timestep with square domain, get the layer-averaged kinetic energy

    Input
    p :: 3-D field psi(z (or m),y,x)

    Output
    KE :: 1-D array of length nz with layer-averaged kinetic energy values
    (total KE in the sytem is given by sum(KE*h))
    """
    
    from scipy.fft import fft2
    
    nz = np.shape(p)[0]
    KE = np.zeros(nz)

    N = np.shape(p)[-1]
    k = np.concatenate([np.arange(0,N/2+1),np.arange(-N/2+1,0)])
    KX, KY = np.meshgrid(k,k)
    K2 = KX**2 + KY**2

    for m in range(nz):

        psi_hat = fft2(p[m,:,:])
        KE[m] = 0.5*np.sum(K2*np.abs(psi_hat)**2)/N**2*(2*np.pi/Lx)**2*Lx*Ly
        
    return KE

def get_TE(p,q,h,Lx=1,Ly=1):
    """
    At a fixed timestep, get total energy

    Input
    p :: 3-D field psi(z (or m),y,x) indexed top to bottom
    q :: 3-D field   q(z (or m),y,x) indexed top to bottom
    h :: layer heights indexed top to bottom

    Output
    TE :: float of total energy in the system
    """
    
    TE = np.sum(np.sum(-0.5*p*q,axis=(1,2))*h)*Lx*Ly
        
    return TE

def get_APE(p,h,S,Lx=1,Ly=1):
    """
    At a fixed timestep, get the available potential energy

    Input
    p :: 3-D field psi(z (or m),y,x) indexed top to bottom
    h :: (nz) array of layer heights indexed top to bottom
    S :: (nz+1) array of stratification (f^2/N^2(zi)) values at cell interfaces

    Output
    APE :: float of available potential energy in the system
    """
    
    nz = np.shape(p)[0]
    APE = np.zeros(np.shape(p))

    APE[0,:,:] = 2*( - S[1]/(h[0]+h[1])*(p[0,:,:]-p[1,:,:]) )*p[0,:,:]

    for k in range(1,nz-1):
        
        APE[k,:,:] = 2*( S[k  ]/(h[k-1]+h[k  ])*(p[k-1,:,:]-p[k  ,:,:]) 
                       - S[k+1]/(h[k  ]+h[k+1])*(p[k  ,:,:]-p[k+1,:,:]) )*p[k,:,:]

    APE[-1,:,:] = 2*( S[-2]/(h[-2]+h[-1])*(p[-2,:,:]-p[-1,:,:]) )*p[-1,:,:]
    
    return -np.sum(APE)/2.*Lx*Ly

def cospectrum(a,b):
    """
    Takes two physical-space real, square fields A and B and calculates their Fourier co-spectrum (adapted from Ian's matlab version)
    
    Input
    a, b :: two square 2-D fields (e.g. p[time,zlevel,:,:])
    
    Output
    cs :: Fourier cospectrum, collapsed into 1-D (i.e. cospectrum over horizontal wavenumber magnitude |k|)
    """
    
    from scipy.fft import fft2, fftfreq
    
    a_hat = fft2(a); b_hat = fft2(b);
    N = np.shape(a)[0]
    
    k = N*fftfreq(N) #np.concatenate([np.arange(0,N/2),np.arange(-N/2,0)])
    KX, KY = np.meshgrid(k,k)
    
    cs = 1j * np.zeros(int(N/2+1))
    for jj in range(N):
        for ii in range(N):
            k = np.sqrt(KX[ii,jj]**2+KY[ii,jj]**2)
            r = k - np.floor(k)
            if k < N/2:
                kk = int(np.floor(k))
                cs[kk]   = cs[kk]  +(1-r)*np.vdot(a_hat[ii,jj],b_hat[ii,jj])
                cs[kk+1] = cs[kk+1]+   r *np.vdot(a_hat[ii,jj],b_hat[ii,jj])
                
    return np.real(cs)/(N**4)

def cospectrum3D(amode,bmode):
    """
    Takes two physical-space real, square fields A and B and calculates their Fourier co-spectrum (adapted from Ian's matlab version)

    Input
    amode, bmode :: two 3-D fields with the vertical (axis 0) projected onto modes (saved as e.g. p_mode.*.dat)

    Output
    cs :: Fourier cospectrum, collapsed into 2-D (vertical m, horizontal |k|)
    """
        
    nmodes = np.shape(amode)[0]
    N = np.shape(amode)[1]
    
    cs = np.zeros((nmodes,int(N/2+1)))
    
    for m in range(nmodes):
        cs[m,:] = cospectrum(amode[m,:,:],bmode[m,:,:])
        
    return cs

def get_TE_spectra(pmode,qmode):
    """
    At a fixed timestep, compute the spectral (vertical m, horizontal |k|) decomposition of energy

    Input
    pmode :: 3-D field psi(mode m,y,x)
    qmode :: 3-D field q(mode m,y,x)

    Output
    2-D array of (nx) x (nx/2+1) of total energy spectral coefficients
    """
    
    return cospectrum3D(-qmode,pmode)

def get_enstrophy_spectra(qmode):
    """
    At a fixed timestep, compute the spectral (vertical m, horizontal |k|) decomposition of enstrophy

    Input
    qmode :: 3-D field q(mode m,y,x)

    Output
    2-D array of (nx) x (nx/2+1) of enstrophy spectral coefficients
    """
            
    return cospectrum3D(qmode,qmode)

## use function to define (sparse) linear operator
from scipy.sparse.linalg import LinearOperator

def Afun(x,M,K):
    
    """
    Define sparse linear operator for the 2D Laplacian with equal unit cells

    Input
    x :: flattened 2D field
    M, K :: 2D dimensions of the field

    Output
    Ax :: flattened result of applying the Laplacian to the input x
    """
    
    
    ## symmetric, positive definite version of the Laplacian. Solve Ax=-b
        
    X = np.reshape(x,(M,K))
    Ax = np.zeros(np.shape(X))

    #### interior ####
    Ax[1:-1,1:-1] = -X[:-2,1:-1] - X[2:,1:-1] - X[1:-1,:-2] - X[1:-1,2:] + 4*X[1:-1,1:-1]

    #### edge boundaries ####
    Ax[ 0,1:-1] = -X[ 1,1:-1] - X[ 0,:-2] - X[ 0,2:] + 3*X[ 0,1:-1]
    Ax[-1,1:-1] = -X[-2,1:-1] - X[-1,:-2] - X[-1,2:] + 3*X[-1,1:-1]
    Ax[1:-1, 0] = -X[:-2, 0] - X[2:, 0] - X[1:-1, 1] + 3*X[1:-1, 0]
    Ax[1:-1,-1] = -X[:-2,-1] - X[2:,-1] - X[1:-1,-2] + 3*X[1:-1,-1]

    #### corner boundaries ####
    Ax[ 0, 0] = -X[0,1] - X[1,0] + 2*X[0,0]
    Ax[ 0,-1] = -X[0,-2] - X[1,-1] + 2*X[0,-1]
    Ax[-1, 0] = -X[-2,0] - X[-1,1] + 2*X[-1,0]
    Ax[-1,-1] = -X[-1,-2] - X[-2,-1] + 2*X[-1,-1]
    
    return Ax.flatten()
