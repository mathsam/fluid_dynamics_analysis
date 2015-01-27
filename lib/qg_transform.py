import numpy as np
import scipy.fftpack as fftpack
import nc_tools

"""
deals with output from qg_model using NetCDF output
assume the field has shape as

    psi(time_step (optional), real_and_imag, ky, kx, z)
"""

def get_vorticity(psik):
    """
    Calculates spectral relative
    vorticity 'zetak' and vortex stretching 'etak' from spectral
    streamfunction field psik
    @param psik stream function in spectral space
                returned from real2complex, complex numpy array
    @return zetak relative vorticity in spectral space
    """
    kmax = psik.shape[-3] - 1
    kx_, ky_ = np.meshgrid(range(-kmax, kmax+1), range(0, kmax+1))
    k2 = kx_**2 + ky_**2
    k2.shape = (1,)*(psik.ndim-3) + psik.shape[-3:-1] + (1,)
    zetak = -k2*psik
    return zetak
    

def real2complex(rfield):
    """
    convert raw qg_model output to complex numpy array
    suppose input has shape
        psi(time_step (optional), real_and_imag, ky, kx, z)
    """
    return rfield[...,0,:,:,:]+1j*rfield[...,1,:,:,:]

def fullspec(hfield):
    """
    Assumes 'hfield' to contain upper-half plane of spectral field, 
    and specifies lower half plane by conjugate 
    symmetry (since physical field is assumed real-valued).  'hfield'
    should have shape (...,kmax+1,2*kmax+1,nz), kmax = 2^n-1,
    hence physical resolution will be 2^(n+1) x 2^(n+1) x nz.  
    NOTE:  The bottom row of the input field corresponds to ky = 0,
    the kx<0 part is NOT assumed a priori to be conjugate-
    symmetric with the kx>0 part.
    """
    if not isinstance(hfield, np.ndarray):
        raise TypeError("input needs to be numpy array")
    if hfield.ndim < 2:
        raise ValueError("array must be at least 2 dimensional")
        
    nky, nkx = hfield.shape[-3:-1]

    if nkx+1 != 2*nky:
        raise ValueError("hfield must have dim (..., kmax+1, 2*kmax+1)")
    hres = nkx + 1
    kmax = nky - 1
    fk = np.zeros(hfield.shape[:-3]+(hres,hres)+(hfield.shape[-1],), 
                  dtype=complex)
    
    fup = np.copy(hfield)
    fup[...,0,kmax-1::-1,:] = fup.conj()[...,0,kmax+1:,:]
    #fup[...,0,kmax,:] = 0. # a littile confused whether should do this
    fdn = np.copy(fup.conj()[..., nky-1:0:-1, nkx-1::-1,:])
    fk[..., nky:, 1:,:] = fup
    fk[...,1:nky, 1:,:] = fdn
    return fk


def spec2grid(sfield):
    """
    Transform one frame of SQG model
    output to a grided (physical) representation.  Assumes 'sfield'
    to be up-half plane, and specifies lower half plane by conjugate
    sym (since physical field is assumed real-valued).  Input field
    should have dimensions  (...,kmax+1,2*kmax+1,nz), where
    kmax=2^n-1, hence physical resolution will be 2^(n+1) x 2^(n+1).
    NOTE: bottom row of the input field corresponds to ky = 0, the
    kx<0 part is NOT assumed a priori to be conjugate- symmetric
    with the kx>0 part.  NOTE: grid2spec(spec2grid(fk)) = fk.
    OPTIONAL: da = true pads input with 0s before transfoming to
    gridspace, for dealiased products.  Default is da = false.
    """
    hres = sfield.shape[-2] + 1
    fk = fullspec(sfield)
    fk = fftpack.ifftshift(fk, axes=(-2,-3))
    return hres*hres*np.real(fftpack.ifft2(fk, axes=(-2,-3)))

def energy_spec(psic):
    """
    Given stream function at one layer of one time step, returns its energy 
    spectrum
    @param psic complex numpy array
                returned from real2complex, has shape (time (optional), ky, kx)
    @return (wavenumber, Ek, EKEk) 1d numpy array
    """
    if not hasattr(energy_spec, 'precal') \
       or (energy_spec.nky, energy_spec.nkx) != psic.shape:
        nky, nkx = psic.shape[-2:]
        energy_spec.nkx   = psic.shape[-1]
        energy_spec.nky   = psic.shape[-2]
        energy_spec.ksqd_ = np.zeros((nky, nkx), dtype=float)
        
        kmax = energy_spec.nky - 1
        for j in range(0, nky):
            for i in range(0, nkx):
                energy_spec.ksqd_[j,i] = (i-kmax)**2 + j**2
        
        radius_arr = np.floor(np.sqrt(energy_spec.ksqd_)).astype(int)
        energy_spec.radius_mask = []
        for i in range(0, kmax):
            energy_spec.radius_mask.append(radius_arr == i+1)
        #print "precalculation done"
    
    kmax = energy_spec.nky - 1
    indx_kx0 = (energy_spec.nkx-1)/2
    KE2d  = np.zeros((energy_spec.nky, energy_spec.nkx), dtype=float)
    Ek    = np.zeros(kmax,       dtype=float)
    EKEk  = np.zeros(kmax,       dtype=float)

    if psic.ndim == 3:
        for i in range(0, psic.shape[0]):
            KE2d += np.abs(psic[i,...]*psic[i,...].conj())
        KE2d *= energy_spec.ksqd_
        KE2d /= psic.shape[0]
    else:
        KE2d  = energy_spec.ksqd_*np.abs(psic*psic.conj())

    for i in range(0,kmax):
        Ek[i]   = np.sum(KE2d[energy_spec.radius_mask[i]])
    
    KE2d[:,indx_kx0] = 0.

    for i in range(0,kmax):
        EKEk[i]   = np.sum(KE2d[energy_spec.radius_mask[i]])

    return np.arange(1,kmax+1), Ek, EKEk
    
def _barotropic_Ek_ncchain(psi):
    """
    process NetCDFChain object with lots of files
    @param psi has shape (time, real_and_imag, ky, kx, z)
    """
    if not isinstance(psi, nc_tools.NetCDFChain):
        raise TypeError("Not NetCDFChain object")
    nky, nkx = psi.shape[-3:-1]
    ksqd_    = np.zeros((nky, nkx), dtype=float)
    kmax     = nky - 1
    indx_kx0 = (nkx-1)/2
    KE2d  = np.zeros((nky, nkx), dtype=float)
    Ek    = np.zeros(kmax,       dtype=float)
    EKEk  = np.zeros(kmax,       dtype=float)
    
    for j in range(0, nky):
        for i in range(0, nkx):
            ksqd_[j,i] = (i-kmax)**2 + j**2
        
    radius_arr = np.floor(np.sqrt(ksqd_)).astype(int)

    for i, dt in enumerate(psi.time_steps_each_file):
        #reading each file all together; if read in each time step each time
        #reading files will take too much time
        psi_seg = psi[psi.time_steps_before[i]:psi.time_steps_before[i]+dt]
        if psi_seg.ndim != 5:
            raise TypeError("file does not have correct number of dimensions")
        psi_seg = np.mean(psi_seg, psi_seg.ndim-1) #barotropic field
        KE2d += np.mean(np.sum(psi_seg*psi_seg, 1), 0)
        
    KE2d *= ksqd_
    KE2d /= len(psi.sorted_files)

    for i in range(0,kmax):
        Ek[i]   = np.sum(KE2d[radius_arr == i+1])
    
    KE2d[:,indx_kx0] = 0.

    for i in range(0,kmax):
        EKEk[i]   = np.sum(KE2d[radius_arr == i+1])

    return np.arange(1,kmax+1), Ek, EKEk
    
    
def barotropic_Ek(psic):
    """
    barotropic energy spectrum for multiple times and multiple layers
    @param psic complex numpy array
                returned from real2complex, has shape 
                (time (optional), ky, kx, z)
            OR  NetCDFChain object
                has shape (time, real_and_imag, ky, kx, z)
    @return (wavenumber, Ek, EKEk) 1d numpy array
    """
    if isinstance(psic, np.ndarray):
        barotropic_psic = np.mean(psic, psic.ndim - 1)
        return energy_spec(barotropic_psic)
    elif isinstance(psic, nc_tools.NetCDFChain):
        return _barotropic_Ek_ncchain(psic)
    else:
        raise TypeError("Input type not supported")