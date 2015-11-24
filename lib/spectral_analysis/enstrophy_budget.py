import qg_transform as qgt
import numpy as np

class EnstrophyTotal(object):
    """barotropic enstrophy budget of the total flow
    Return time and zonally averaged
    """
    def __init__(self, psic):
        """
        Args:
            psic: stream function complex spectrum, expect the shape to be
                (time, ky, kx, z). nz = 2
        """
        self._bt_psis = 0.5*(psic[...,0] + psic[...,1])
        self._bc_psis = 0.5*(psic[...,0] - psic[...,1])
        self._vorg = qgt.spec2grid(qgt.get_vorticity(self._bt_psis))
    
    def advection(self):
        """Evaluate -vor * J_kl(psi, vor))
        """
        bt_vors = qgt.laplacian(self._bt_psis, 1)
        Jg_psi_vor = qgt.jacobian(self._bt_psis, bt_vors)
        return -np.mean(np.mean(self._vorg*Jg_psi_vor, -1), 0)
    
    def enstrophy_generation(self):
        """ -vor * J_kl(tau, \nabla tau))
        """
        Jg_tau_bcvor = qgt.jacobian(self._bc_psis, qgt.laplacian(self._bc_psis))
        return -np.mean(np.mean(self._vorg*Jg_tau_bcvor, -1), 0)
    
    def enstrophy_generation_shear(self, U = 1):
        """ -vor * \nabla dtau/dx)
        
        Args:
            U: mean shear parameter
        """
        nabla_dx_taug = qgt.spec2grid(qgt.laplacian(qgt.partial_x(self._bc_psis)))
        return -U * np.mean(np.mean(self._vorg*nabla_dx_taug, -1), 0)
    
class PETotal(object):
    """Potential energy budget for the total flow
    """
    def __init__(self, psic, kd):
        self._kd = kd
        self._bt_psis = 0.5*(psic[...,0] + psic[...,1])
        self._bc_psis = 0.5*(psic[...,0] - psic[...,1])
        self._taug = qgt.spec2grid(qgt.get_vorticity(self._bc_psis))
    
    def bc2bt(self):
        Jg_psi_bcvor = qgt.jacobian(self._bt_psis, qgt.laplacian(self._bc_psis))
        return np.mean(np.mean(self._taug*Jg_psi_bcvor, -1), 0)
    
    def pe_adv(self):
        v = qgt.get_velocities(self._bt_psis)[1]
        vg = qgt.spec2grid(v)
        dy_taug = qgt.spec2grid(qgt.partial_y(self._bc_psis))
        return -self._kd**2 * np.mean(np.mean(self._taug*vg*dy_taug, -1), 0)