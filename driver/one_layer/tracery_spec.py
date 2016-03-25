import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np

filename_prefix = 'Dec12_kf16_drag64e-4'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
psi = qg_transform.real2complex(
    nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                         '%s_seg' %filename_prefix,'psi', last_n_files=1))
tracery = qg_transform.real2complex(
    nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                         '%s_seg' %filename_prefix,'tracer_y', last_n_files=1))                      
v = qg_transform.get_velocities(psi)[1]
k, tracer_flux_k = qg_transform.prod_spectrum(v, tracery)
k, v2_k          = qg_transform.prod_spectrum(v, v)
k, tracer2_k     = qg_transform.prod_spectrum(tracery, tracery)
##
plt.semilogx(k, -tracer_flux_k_kf16/np.sqrt(v2_k_kf16*tracer2_k_kf16), label='kf=16, drag=256e-4')
plt.semilogx(k, -tracer_flux_k_kf16_drag64/np.sqrt(v2_k_kf16_drag64*tracer2_k_kf16_drag64), label='kf=16, drag=64e-4')
plt.legend(loc='best')
plt.ylabel('tracer v corr')
plt.show()
##
plt.loglog(k, -tracer_flux_k_kf16, label='kf=16')
plt.loglog(k, -tracer_flux_k_kf32, label='kf=32')
plt.loglog(k, -tracer_flux_k_kf64, label='kf=64')
plt.legend(loc='best')
plt.ylabel('tracer flux')
plt.show()
##
plt.semilogx(k, -tracer_flux_k_kf16/np.sqrt(v2_k_kf16*tracer2_k_kf16), label='kf=16')
plt.semilogx(k, -tracer_flux_k_kf32/np.sqrt(v2_k_kf32*tracer2_k_kf32), label='kf=32')
plt.semilogx(k, -tracer_flux_k_kf64/np.sqrt(v2_k_kf64*tracer2_k_kf64), label='kf=64')
plt.legend(loc='best')
plt.ylabel('tracer v corr')
plt.show()
##
plt.loglog(k, -tracer_flux_k_kf16/np.sqrt(v2_k_kf16), label='kf=16')
plt.loglog(k, -tracer_flux_k_kf32/np.sqrt(v2_k_kf32), label='kf=32')
plt.loglog(k, -tracer_flux_k_kf64/np.sqrt(v2_k_kf64), label='kf=64')
x = np.arange(10, 250).astype(float)
plt.loglog(x, 1e-2/x, label=r'$k^{-1}$')
plt.loglog(x, 1e-1/x**2, label=r'$k^{-2}$')
plt.legend(loc='best')
plt.ylabel('tracer_flux/v')
plt.show()