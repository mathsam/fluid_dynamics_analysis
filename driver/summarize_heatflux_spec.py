import os
import numpy as np
import diag
import matplotlib.pyplot as plt
import matlab_style as ml
analysis_dir = '/home/j1c/analysis/2015/qg_model'
save_dir = '/home/j1c/analysis/2015/qg_model/heat_flux'

exp_list  = ['Jan17_drag_', 'Jan18_c2_drag_', 'Jan18_c2.5_drag_']
cri_list  = ['1.6', '2.0', '2.5']
drag_list = ['1e0', '1e-1', '1e-2', '1e-3', '5e-4', '1e-4']
num_wavenumbers = 511

flux_spec = np.zeros((len(cri_list), len(drag_list), num_wavenumbers))
v2_spec = np.zeros((len(cri_list), len(drag_list), num_wavenumbers))
tau2_spec = np.zeros((len(cri_list), len(drag_list), num_wavenumbers))
v_tau_corr = np.zeros((len(cri_list), len(drag_list), num_wavenumbers))

for i, exp_i in enumerate(exp_list):
    for j, drag_i in enumerate(drag_list):
        exp_name = exp_i + drag_i
        file = os.path.join(analysis_dir, exp_name, 'flux_spec.npz')
        spec_zip = np.load(file)        
        
        flux_spec[i,j,:] = spec_zip['tau_flux'].flatten()
        v2_spec[i,j,:] = spec_zip['v2'].flatten()
        tau2_spec[i,j,:] = spec_zip['tau2'].flatten()
        
        
v_tau_corr = flux_spec/np.sqrt(v2_spec*tau2_spec)
ml.save(os.path.join(save_dir, 'Jan17_18_heatflux.npz'), flux_spec=flux_spec,
                                                    v2_spec=v2_spec,
                                                    tau2_spec=tau2_spec,
                                                    v_tau_corr=v_tau_corr)
## plot correlation
plt.semilogx(v_tau_corr[2,1::2,:].T)
plt.xlabel('Wavenumber')
plt.ylabel(r"corr($\tau'$, $v'$)")
plt.show()
## plot heat flux
plt.loglog(flux_spec[2,1::2,:].T)
plt.xlabel('Wavenumber')
plt.ylabel(r"$\tau'v$")
plt.show()