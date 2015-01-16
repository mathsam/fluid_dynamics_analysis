import numpy as np
import scipy.fftpack as fftpack

def fullspec(hfield):
    """
    Assumes 'hfield' to contain upper-half plane of spectral field, 
    and specifies lower half plane by conjugate 
    symmetry (since physical field is assumed real-valued).  'hfield'
    should have shape (2*kmax+1,kmax+1,:,:), kmax = 2^n-1,
    hence physical resolution will be 2^(n+1) x 2^(n+1) x nz.  
    NOTE:  The bottom row of the input field corresponds to ky = 0,
    the kx<0 part is NOT assumed a priori to be conjugate-
    symmetric with the kx>0 part.
    """
    if not isinstance(hfield, np.ndarray):
        raise TypeError("input needs to be numpy array")
    if hfield.ndim < 2:
        raise ValueError("array must be at least 2 dimensional")
    
    nkx, nky = hfield.shape[0:2]
    if nkx+1 != 2*nky:
        raise ValueError("hfield must have dim (1:2*kmax+1,1:kmax+1,...)")
    hres = nkx + 1
    kmax = nky - 1
    fk = np.zeros((hres,hres)+hfield.shape[2:], dtype=complex)
    
    fup = np.copy(hfield)
    fup[kmax-1::-1, 0,...] = fup.conj()[kmax+1:,0,...]
    #fup[kmax,0,...] = 0. # a littile confused whether should do this
    fdn = np.copy(fup.conj()[nkx-1::-1, nky-1:0:-1, ...])
    fk[1:, nky:, ...] = fup
    fk[1:, 1:nky,...] = fdn
    return fk


def spec2grid(sfield):
    """
    Transform one frame of SQG model
    output to a grided (physical) representation.  Assumes 'sfield'
    to be up-half plane, and specifies lower half plane by conjugate
    sym (since physical field is assumed real-valued).  Input field
    should have dimensions (1:2*kmax+1,1:kmax+1,:,:), where
    kmax=2^n-1, hence physical resolution will be 2^(n+1) x 2^(n+1).
    NOTE: bottom row of the input field corresponds to ky = 0, the
    kx<0 part is NOT assumed a priori to be conjugate- symmetric
    with the kx>0 part.  NOTE: grid2spec(spec2grid(fk)) = fk.
    OPTIONAL: da = true pads input with 0s before transfoming to
    gridspace, for dealiased products.  Default is da = false.
    """
    hres = sfield.shape[0] + 1
    fk = fullspec(sfield)
    fk = fftpack.ifftshift(fk, axes=(0,1))
    return hres*hres*np.real(fftpack.ifft2(fk, axes=(0,1)))
