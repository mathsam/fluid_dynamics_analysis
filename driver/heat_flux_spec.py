import qg_transform
import nc_tools
import numpy as np
import diag
import matplotlib.pyplot as plt

# Sc = 1.6, F = 1583.14, beta = 3957.86
# Sc = 2.0, F = 1583.14, beta = 3166.29
# Sc = 2.5, F = 1583.14, beta = 2533.03
F    = 1583.14
beta = 2533.03

filename_prefix = 'Jan17_drag_1e-4'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
psi = nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                      '%s_seg' %filename_prefix,'psi', last_n_files=1)[:]
psi = qg_transform.real2complex(psi)

## evaluate the spectrum
import matlab_style as ml
v = qg_transform.get_velocities(psi)[1]
bt_v = np.mean(v, -1)[..., np.newaxis]
tau = 0.5*(psi[...,0] - psi[...,1])[..., np.newaxis]
tau = qg_transform.filter(tau, None, None,True)

k, v2 = qg_transform.prod_spectrum(bt_v, bt_v)
k, tau2 = qg_transform.prod_spectrum(tau, tau)

k, tau_flux = qg_transform.prod_spectrum(bt_v, tau)

ml.save(save_dir + 'flux_spec', k=k, v2=v2, tau2=tau2, tau_flux=tau_flux)
## plot spectrum of v'2, tau^2, heat flux
most_energetic_k = np.where(v2 == v2.max())[0][0]
x   = np.arange(10,101)
p53  = v2[most_energetic_k]*x**(-5./3.)/(x[most_energetic_k]**(-5./3.))/5.0

plt.loglog(k, v2, label=r"$v'$")
plt.loglog(k, tau2, label=r"$\tau'$")
plt.loglog(k, tau_flux, label='heat flux')
plt.loglog(x, p53, label='$k^{5/3}$')
plt.legend(loc='best')
plt.xlabel('Wavenumber')
plt.savefig(save_dir + 'heat_flux_spec.png')
plt.show()