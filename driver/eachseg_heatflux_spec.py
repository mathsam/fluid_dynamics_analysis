import qg_transform
import nc_tools
import numpy as np
import diag
import matplotlib.pyplot as plt
import matlab_style as ml

filename_prefix = 'Jan17_drag_1e-4'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix

for i in range(1, 181):
    psi = nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                        '%s_seg%d.nc' %(filename_prefix, i),'psi')
    psi = qg_transform.real2complex(psi)
    v = qg_transform.get_velocities(psi)[1]
    bt_v = np.mean(v, -1)[..., np.newaxis]
    tau = 0.5*(psi[...,0] - psi[...,1])[..., np.newaxis]
    tau = qg_transform.filter(tau, None, None,True)
    
    k, v2 = qg_transform.prod_spectrum(bt_v, bt_v)
    k, tau2 = qg_transform.prod_spectrum(tau, tau)
    
    k, tau_flux = qg_transform.prod_spectrum(bt_v, tau)
    
    ml.save(save_dir + 'flux_spec_seg%d' %i, k=k, v2=v2, tau2=tau2, tau_flux=tau_flux)
    
    most_energetic_k = np.where(v2 == v2.max())[0][0]
    x   = np.arange(10,101)
    p53  = v2[most_energetic_k]*x**(-5./3.)/(x[most_energetic_k]**(-5./3.))/5.0
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog(k, v2, label=r"$v'$")
    ax.loglog(k, tau2, label=r"$\tau'$")
    ax.loglog(k, tau_flux, label='heat flux')
    ax.loglog(x, p53, label='$k^{5/3}$')
    ax.legend(loc='best')
    ax.set_xlabel('Wavenumber')
    fig.savefig(save_dir + 'heatflux_spec_seg%d.png' %i)