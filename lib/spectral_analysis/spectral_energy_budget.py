"""Spectral KE analysis for 2-layer QG model
See formulation in documentation"""
import qg_transform as qgt

class BarotropicKE(object):
    """barotropic total flow
    
    Methods return (k, dKEk/dt of various terms)
    """
    def __init__(self, psic):
        """
        Args:
            psic: stream function complex spectrum, expect the shape to be
                (time, ky, kx, z). nz = 2
        """
        self._bt_psis = 0.5*(psic[...,0] + psic[...,1])
        self._bc_psis = 0.5*(psic[...,0] - psic[...,1])
    
    def nonlinear_transfer(self):
        """Evaluate Re(psi_kl' * J_kl(psi, vor)), ' denote complex conj
        Term I in formulation
        """
        bt_vors = qgt.laplacian(self._bt_psis, 1)
        Jg_psi_vor = qgt.jacobian(self._bt_psis, bt_vors)
        Js_psi_vor = qgt.grid2spec(Jg_psi_vor)
        return qgt.prod_spectrum(self._bt_psis, Js_psi_vor)
    
    def ke_generation(self):
        """Re(psi_kl' * J_kl(tau, \nabla tau))
        Term II in formulation
        """
        Js_tau_bcvor = qgt.grid2spec(qgt.jacobian(self._bc_psis,
            qgt.laplacian(self._bc_psis)))
        return qgt.prod_spectrum(self._bt_psis, Js_tau_bcvor)
    
    def ke_generation_shear(self, U = 1):
        """U*Re(psi_kl' * \nabla dtau/dx)
        Term III in formulation
        
        Args:
            U: mean shear parameter
        """
        nabla_dx_tau = U * qgt.laplacian(qgt.partial_x(self._bc_psis))
        return qgt.prod_spectrum(self._bt_psis, nabla_dx_tau)
    
    def friction_diss(self, kf = 1):
        """kf*Re(psi_kl' * \nabla(psi-tau)/2)
        
        Args:
            kf: frictional coefficient of the bottom layer
        """
        vors1 = kf * 0.5*qgt.laplacian(self._bt_psis - self._bc_psis)
        return qgt.prod_spectrum(self._bt_psis, vors1)
    
    def hyper_diss(self, order=4, mu=1):
        """mu * Re(psi_kl' * \nabla^4(\nabla psi))
        
        Args:
            order: hyperviscosity order
            mu: hyperviscosity coeff
        """
        del_vors = mu * qgt.laplacian(self._bt_psis, 1+order)
        return qgt.prod_spectrum(self._bt_psis, del_vors)

class BarotropicEKE(object):
    """barotropic eddy flow spectral budget
    
    Methods return (k, dEKEk/dt of various terms)
    """
    def __init__(self, psic):
        """
        Args:
            psic: stream function complex spectrum, expect the shape to be
                (time, ky, kx, z). nz = 2
        """
        self._bt_psis = 0.5*(psic[...,0] + psic[...,1])
        self._bc_psis = 0.5*(psic[...,0] - psic[...,1])
    
    @staticmethod
    def _jacob_eddy_mean(psis, zetas):
        """eddy mean interaction part: 
        v'*partial_y(zeta_mean) + u_bar*partial_x(zeta')
        
        Args:
            psis, zetas: spectral fields with shape (time, ky, kx, z)
        
        Return:
            spectrum of the eddy mean interaction field
        """
        zetas_mean = qgt.get_zonalmean(zetas)
        dy_zetag_mean = qgt.spec2grid(qgt.partial_y(zetas_mean))
        del zetas_mean
        dx_zetag_eddy = qgt.spec2grid(qgt.partial_x(zetas))
        us, vs = qgt.get_velocities(psis)
        us = qgt.get_zonalmean(us)
        vg = qgt.spec2grid(vs)
        del vs
        ug = qgt.spec2grid(us)
        del us
        return qgt.grid2spec(vg*dy_zetag_mean + ug*dx_zetag_eddy)
    
    @staticmethod
    def _jacob_eddy_eddy(psis, zetas):
        """eddy-eddy interaction:
        J(psi', zeta')
        
        Args:
            psis, zetas: spectral fields
        Return:
            spectrum of eddy eddy interaction field
        """
        psis_eddy  = qgt.get_eddy(psis)
        zetas_eddy = qgt.get_eddy(zetas)
        return qgt.grid2spec(qgt.jacobian(psis_eddy, zetas))
    
    def eddy_eddy_transfer(self):
        return qgt.prod_spectrum(qgt.get_eddy(self._bt_psis), 
            BarotropicEKE._jacob_eddy_eddy(self._bt_psis, qgt.laplacian(self._bt_psis)))
    
    def eddy_mean_transfer(self):
        return qgt.prod_spectrum(self._bt_psis,
            BarotropicEKE._jacob_eddy_mean(self._bt_psis, qgt.laplacian(self._bt_psis)))
    
    def eddy_eddy_gen(self):
        return qgt.prod_spectrum(qgt.get_eddy(self._bt_psis),
            BarotropicEKE._jacob_eddy_eddy(self._bc_psis, qgt.laplacian(self._bc_psis)))
    
    def eddy_mean_gen(self):
        return qgt.prod_spectrum(self._bt_psis,
            BarotropicEKE._jacob_eddy_mean(self._bc_psis, qgt.laplacian(self._bc_psis)))