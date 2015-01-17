import numpy as np
import scipy.fftpack as fftpack

def get_vorticity(psik):
    """
    Calculates spectral relative
    vorticity 'zetak' and vortex stretching 'etak' from spectral
    streamfunction field psik
    @param psik stream function in spectral space
    @return zetak relative vorticity in spectral space
    """
    kmax = psik.shape[-2] - 1
    kx_, ky_ = np.meshgrid(range(-kmax, kmax+1), range(0, kmax+1))
    k2 = kx_**2 + ky_**2
    k2.shape = (1,)*(psik.ndim-2) + psik.shape[-2:]
    zetak = -k2*psik
    return zetak
    

def real2complex(rfield):
    """
    convert qg_model output to complex numpy array
    suppose the last dimension is real_and_imag
    """
    return rfield[:,0]+1j*rfield[:,1]

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
