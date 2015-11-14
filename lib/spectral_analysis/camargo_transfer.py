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
    
    def is_valid_mode(self, kx, ky):
        if abs(kx) > self._kxymax or abs(ky) > self._kxymax:
            return False
        return True
    
    def get_mode(self, kx, ky):
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
    def __init__(self, spec, pre_cal=True, K_wall=50):
        """
        Args:
            spec: StandardSpectrum object or numpy array (not std spec)
            pre_cal: bool|numpy array, whether to precalculate all the transfers within
                      wavenumber shell `K_wall`
                     if type is numpy array, just assign this array to be _transM
            K_wall: integer
        """
        if isinstance(spec, StandardSpectrum):
            self._stdspec = spec
        else:
            self._stdspec = StandardSpectrum(spec)
        
        if isinstance(pre_cal, bool):
            if pre_cal:
                num_ks = 2*K_wall + 1
                self._transM = np.zeros((num_ks, num_ks, num_ks, num_ks)) #(kx1, ky1, kx, ky)
                kx1_set = range(-K_wall, K_wall+1)
                ky1_set = range(-K_wall, K_wall+1)
                kx_set  = range(-K_wall, K_wall+1)
                ky_set  = range(-K_wall, K_wall+1)
                for i1, kx1 in enumerate(kx1_set):
                    for j1, ky1 in enumerate(ky1_set):
                        print "calculating kx1=%d, ky1=%d" %(kx1, ky1)
                        for i, kx in enumerate(kx_set):
                            for j, ky in enumerate(ky_set):
                                self._transM[i1,j1,i,j] = CamargoTrans._cal_transfer(
                                    self._stdspec, (kx1, ky1), (kx, ky))
                self._if_pre_cal = True
                self._K_wall = K_wall
        elif isinstance(pre_cal, np.ndarray):
            num_kxs, num_kys = pre_cal.shape[:2]
            if num_kxs == num_kys:
                self._transM = pre_cal.copy()
                self._K_wall = (num_kxs-1)/2
            else:
                raise IndexError('pre_cal array has wrong shape')
        else:
            raise TypeError('pre_cal has wrong type')
    
    @staticmethod
    def _cal_transfer(stdspec, K1, K):
        """
        Args:
            stdspec: StandardSpectrum object
            K1: tuple for (k1x, k1y)
            K: tuple for (kx, ky)
        
        Return:
            a float
        """
        k1x, k1y = K1
        kx, ky   = K
        k2x = kx - k1x
        k2y = ky - k1y
        
        if not stdspec.is_valid_mode(k2x, k2y):
            raise IndexError('interaction not exist')
        
        psi_K1 = stdspec.get_mode(k1x, k1y)
        psi_K2 = stdspec.get_mode(k2x, k2y)
        psi_K  = stdspec.get_mode(kx, ky)
        
        return -2.*(ky*k1x - k1y*kx)*(k2x**2+k2y**2)*np.mean(
            np.real(np.conj(psi_K) * psi_K2 * psi_K1) )
        
    
    def get_transfer(self, K1, K):
        """
        Args:
            K1: tuple for (k1x, k1y)
            K:  tuple for (kx, ky)
        
        Return:
            a float
        """
        if hasattr(self, '_transM'):
            k1x, k1y = K1
            kx, ky   = K
            return self._transM[k1x+self._K_wall, k1y+self._K_wall,
                                kx+self._K_wall,  ky+self._K_wall]
        return CamargoTrans._cal_transfer(self._stdspec, K1, K)