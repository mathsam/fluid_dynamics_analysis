import sys
sys.path.append('/home/j1c/py/lib')
sys.path.append('/home/j1c/py/lib/spectral_analysis')
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style as ml
import barotropic_spec_diag as bt_sp

filename_prefix = 'Dec12_kf128_drag8e-4'
filename = r'%s_seg150' %filename_prefix

filedir  = '/archive/Junyi.Chai/QG_exp/2015/%s' %filename_prefix
beta = 16*np.pi

print "working on %s" %filename
psic = qg_transform.real2complex(nc_tools.ncread(filedir, filename,'psi'))
tracerc = qg_transform.real2complex(nc_tools.ncread(filedir, filename,'tracer_y'))
vorc = qg_transform.get_vorticity(psic)

diff_vor_t = vorc/beta-tracerc
k, E_vor_t_k = qg_transform.prod_spectrum(diff_vor_t, diff_vor_t)

##
tracerc = qg_transform.filter(tracerc, None, None, True)
vorc = qg_transform.filter(vorc, None, None, True)
k, E_t_vor = qg_transform.prod_spectrum(tracerc, vorc)
k, E_t = qg_transform.prod_spectrum(tracerc, tracerc)
k, E_vor = qg_transform.prod_spectrum(vorc, vorc)
corr_t_vor = E_t_vor/np.sqrt(E_t*E_vor)
##
plt.semilogx(k, corr_t_vor_kf16,label='kf=16')
#plt.semilogy(k, corr_t_vor_kf32,label='kf=32')
plt.semilogx(k, corr_t_vor_kf64,label='kf=64')
plt.semilogx(k, corr_t_vor_kf128,label='kf=128')
plt.plot([16, 16], [1, 0.2], '--k')
plt.plot([64, 64], [1, 0.2], '--k')
plt.plot([128, 128], [1, 0.2], '--k')
plt.text(16, 0.1, r'$k_f=16$')
plt.text(64, 0.1, r'$k_f=64$')
plt.text(128, 0.1, r'$k_f=128$')

plt.plot([24, 24], [1, 0.4], '--g')
plt.text(24, 0.8, r'$k_\beta$')
plt.ylabel(r"$corr(\zeta', S')$")
plt.xlabel('Wavenumber')
plt.legend(loc='best')
plt.show()