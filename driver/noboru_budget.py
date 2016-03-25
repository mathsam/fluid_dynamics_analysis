import sys
sys.path.append('/home/j1c/py/lib')
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style

filename_prefix = 'Nov5_Sc1.6_drag5e-2'
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
filename = r'%s_seg[0-9]+' %filename_prefix
psic = qg_transform.real2complex(
    nc_tools.NetCDFChain(filedir, filename,'psi',last_n_files=1))

params = {'beta': 3957.85873603000, 
          'F'   : 1583.14349441000, 
          'bot_drag' : 3.978873577300000E}#3.978873577300000E-002} 

fwd = TwoLayerWaveActivity(psic, params)
fwd.cal_barotropic()
#fwd.cal_upperlayer()
#fwd.cal_lowerlayer()
vorflux = vorflux_1d(np.mean(psic, -1))
# westerly jet
bt_psi = np.mean(psic, -1)
u, _ = qg_transform.get_velocities(bt_psi)

u_ave = np.mean(u, 0)[np.newaxis,...]
ug_ave = np.squeeze(qg_transform.spec2grid(u_ave))
ug1d = np.mean(ug_ave, -1)

lat_grid = np.arange(0, 1024)
#westerly_idx = np.logical_and(np.logical_and(ug1d>0, lat_grid>300),
#    lat_grid<450)
westerly_idx = ug1d>0
## save calculations
hypvis_btpv = fwd.hypvis_btpv.copy()
mU_vor_dtaudx = fwd.mU_vor_dtaudx.copy()
mJ_tau_vortau = fwd.mJ_tau_vortau.copy()
drag_btpv = fwd.drag_btpv.copy()
lang_btpv = fwd.lang_btpv.copy()

save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
import os
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)
matlab_style.save(save_dir + 'noboru_budget', 
    hypvis_btpv = hypvis_btpv, mU_vor_dtaudx = mU_vor_dtaudx,
    mJ_tau_vortau = mJ_tau_vortau, drag_btpv = drag_btpv,
    lang_btpv = lang_btpv,
    vorflux = vorflux, ug1d = ug1d)
    
##

plt.plot(hypvis_btpv, label='hypvis')
plt.plot(mU_vor_dtaudx+mJ_tau_vortau, label='enstrophy forcing')
plt.plot(drag_btpv, label='drag')
plt.plot(hypvis_btpv+mU_vor_dtaudx+mJ_tau_vortau, label='hyp+forc')
plt.plot(hypvis_btpv+mU_vor_dtaudx+mJ_tau_vortau+drag_btpv, label='total')
plt.plot(vorflux, label='vorflux')
plt.legend(loc='best')
plt.show()

## save noboru_budget
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
import os
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

save_name = 'noboru_budget'
plt.savefig(save_dir + save_name)
##
plt.scatter(dPVdy[westerly_idx], forc[westerly_idx], label='forc')
plt.scatter(dPVdy[westerly_idx], -hypvis_btpv[westerly_idx], c='r', label='diffusive flux')
plt.scatter(dPVdy[westerly_idx], forc[westerly_idx]+fwd.hypvis_btpv[westerly_idx], c='g', label='forc+diff_flux')
plt.legend(loc='best')
plt.show()
#plt.savefig(save_dir + 'dPVdy_vs_forc_diffflux')
## regression: dPVdy vs forc
import scipy.stats as stats
from scipy.optimize import curve_fit
dPVdy = dfdy(lang_btpv)
forc  = fwd.mU_vor_dtaudx+fwd.mJ_tau_vortau
x = dPVdy[westerly_idx]
y = forc[westerly_idx]
F0 = curve_fit(lambda x, K: np.log(K/x), x, np.log(y), (x*y).mean())[0][0]
k, b = stats.linregress(np.log10(x), np.log10(y))[:2]
plt.scatter(x, forc[westerly_idx], label='dPVdy vs forc')
x1 = np.linspace(np.min(x), np.max(x), 100)
plt.plot(x1, F0/x1, label='%.1e/x' %F0)
plt.plot(x1, x1**(k) * 10**b, label='y=%.1e*x^%1.e' %(10**b, k))
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.xlabel('dPV/dy')
plt.ylabel('forc')
plt.legend(loc='best')
save_name = 'dPVdy_vs_forc_loglog'
plt.show()
## regression: loglog diffusive flux vs dPVdy
from scipy.optimize import curve_fit
x = dPVdy[westerly_idx].copy()
y = -hypvis_btpv[westerly_idx].copy()/x
#K0, mu_d_beta2d = curve_fit(lambda x, K0, mu_d_beta2: K0/(1+mu_d_beta2*x**2), x, y,
#    [30./25000., 6./2203**2])[0]
K0, mu_d_beta2d = curve_fit(lambda x, K0, mu_d_beta2: x*K0/(1+mu_d_beta2*x**2), x, x*y,
    [30./25000., 6./2203**2])[0]
k, b = stats.linregress(np.log10(x), np.log10(x*y))[:2]
plt.scatter(x, x*y, label='dPVdy vs diffusive flux')
x1 = np.linspace(np.min(x), np.max(x), 100)
plt.plot(x1, x1*K0/(1+mu_d_beta2d*x1**2), label=r'$K_0$=%.1e, $\mu\beta^{-2}$=%.1e' %(K0, mu_d_beta2d))
plt.plot(x1, x1**(k) * 10**b, label='y=%.1e*x^%1.e' %(10**b, k))
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.legend(loc='best')
plt.xlabel('dPV/dy')
plt.ylabel('diffusive PV flux')
#plt.show()
save_name = 'dPVdy_vs_diffflux_loglog'

## regression: diffusive flux vs dPVdy
from scipy.optimize import curve_fit
x = dPVdy[westerly_idx].copy()
y = -fwd.hypvis_btpv[westerly_idx].copy()/x
#K0, mu_d_beta2d = curve_fit(lambda x, K0, mu_d_beta2: K0/(1+mu_d_beta2*x**2), x, y,
#    [30./25000., 6./2203**2])[0]
K0, mu_d_beta2d = curve_fit(lambda x, K0, mu_d_beta2: x*K0/(1+mu_d_beta2*x**2), x, x*y,
    [30./25000., 6./2203**2])[0]

plt.scatter(x, x*y, label='dPVdy vs diffusive flux')
x1 = np.linspace(np.min(x), np.max(x), 100)
plt.plot(x1, x1*K0/(1+mu_d_beta2d*x1**2), label=r'$K_0$=%.1e, $\mu\beta^{-2}$=%.1e' %(K0, mu_d_beta2d))
plt.legend(loc='best')
plt.xlabel('dPV/dy')
plt.ylabel('diffusive PV flux')
plt.show()

## save figure
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
import os
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

#save_name = 'dPVdy_vs_diffflux_loglog'
plt.savefig(save_dir + save_name)
##
vg = qg_transform.spec2grid(qg_transform.get_velocities(psic)[1])
pvflux = np.mean(np.mean(vg*finite_wave_diag._pvg, 0), 1)