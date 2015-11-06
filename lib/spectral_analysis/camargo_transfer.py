import numpy as np

class StandardSpectrum(object):
    def __init__(self, spec):
        """wrap around spectrum field `spec` organized as in Shafer's qg model
        in order to use it in a more standard fashion
        
        Args:
            spec: complex numpy array with shape (time(optional), ky, kx)
        """
        self._spec = spec
        self._kxymax = spec.shape[-2]-1
    
    def is_valid_mode(kx, ky):
        if abs(kx) > self._kxymax or abs(ky) > self._kxymax:
            return False
        return True
    
    def get_mode(kx, ky):
        """get complex component for mode (kx, ky)
        
        Args:
            kx, ky: integers
        """
        if abs(kx) > self._kxymax or abs(ky) > self._kxymax:
            raise IndexError('mode (%d, %d) does not exist' %(kx, ky) )
        
        take_conj = False
        if ky < 0:
            ky = -ky
            take_conj = True
        elif ky == 0 and kx < 0:
            kx = -kx
            take_conj = True
            
        if take_conj:
            return np.conj(self._spec[..., ky, kx+self._kxymax])
        else:
            return self._spec[..., ky, kx+self._kxymax]
            
class CamargoTrans(object):
    """energy transfer from mode (k1x, k1y) to (kx, ky)
    
    See Manz, Ramisch and Stroth, 2009, Physical mechanism behind zonal-flow
    generation in drift-wave turbulence, Phys. Rev. Lett.
    """
    def __init__(self, spec):
        """
        Args:
            spec: StandardSpectrum object or numpy array (not std spec)
        """
        if isinstance(spec, StandardSpectrum):
            self._stdspec = spec
        else:
            self._stdspec = StandardSpectrum(spec)
    
    def get_transfer(K1, K):
        """
        Args:
            K1: tuple for (k1x, k1y)
            K:  tuple for (kx, ky)
        
        Return:
            a float
        """
        k1x, k1y = K1
        kx, ky   = K
        k2x = kx - k1x
        k2y = ky - k1y
        
        if not self._stdspec.is_valid_mode(k2x, k2y):
            raise IndexError('interaction not exist')
        
        psi_K1 = self._stdspec.get_mode(k1x, k1y)
        psi_K2 = self._stdspec.get_mode(k2x, k2y)
        psi_K  = self._stdspec.get_mode(kx, ky)
        
        return -2.*(ky*k1x - k1y*kx)*(k2x**2+k2y**2)*np.mean(
            np.real(np.conj(psi_K) * psi_K2 * psi_K1)) )