import sys
sys.path.append('/home/j1c/py/lib')
sys.path.append('/home/j1c/py/lib/noboru')
from finite_amp_diag import OneLayerWaveActivity, vorflux_1d, dfdy
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style

filename_prefix = 'Dec12_kf16_drag8e-4'
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
filename = r'%s_seg[0-9]+' %filename_prefix
psic_bt  = qg_transform.real2complex(
    nc_tools.NetCDFChain(filedir, filename,'psi',last_n_files=1))[...,0]
rm_forcing = qg_transform.real2complex(
    nc_tools.NetCDFChain(filedir, filename,'rm_forcing',last_n_files=1))

params = {'beta': 50.2654824574000, 
          'bot_drag' : 256e-4}

fwd = OneLayerWaveActivity(psic_bt, params)
fwd.cal_barotropic(rm_forcing=rm_forcing, num_equlats=1024, zoom=1, border_y=256)

u, v = qg_transform.get_velocities(psic_bt)

u_ave = np.mean(u, 0)[np.newaxis,...]
ug_ave = np.squeeze(qg_transform.spec2grid(u_ave))
ug1d = np.mean(ug_ave, -1)

lat_grid = np.arange(0, 1024)
#westerly_idx = np.logical_and(np.logical_and(ug1d>0, lat_grid>300),
#    lat_grid<450)
westerly_idx = ug1d>0

vorflux = vorflux_1d(psic_bt)
##
tracer_y = qg_transform.real2complex(
    nc_tools.NetCDFChain(filedir, filename,'tracer_y',last_n_files=1))[...,0]
tracer_yg = qg_transform.spec2grid(tracer_y)
lang_ty, A_ty = time_mean_wave_activity_op(tracer_yg, 1024)
beta2d = params['beta']*np.linspace(-np.pi, np.pi, 1024)[:,np.newaxis]
##
hypvis_btpv = fwd.hypvis_btpv.copy()
drag_btpv = fwd.drag_btpv.copy()
lang_btpv = fwd.lang_btpv.copy()
rm_forc   = fwd.rm_forc.copy()
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
import os
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

matlab_style.save(save_dir + 'noboru_budget_zoom1', 
    hypvis_btpv = hypvis_btpv, drag_btpv = drag_btpv,
    lang_btpv = lang_btpv,
    ug1d = ug1d,
    rm_forc=rm_forc,
    vorflux=vorflux)
##
#plt.plot(ug1d*hypvis_btpv.std()/ug1d.std(), '--', label='zonal wind')
plt.plot(hypvis_btpv, label='hypvis')
plt.plot(rm_forc, label='enstrophy forcing')
plt.plot(drag_btpv, label='drag')
plt.plot(hypvis_btpv+rm_forc+drag_btpv, label='total')
#plt.plot(vorflux, label='vorflux')
plt.plot(params['bot_drag']*ug1d, label='bottom_drag')
plt.legend(loc='best')
plt.show()

## fit dPV/dy vs. diffusive flux
import scipy.stats as stats
from scipy.optimize import curve_fit
dPVdy = dfdy(lang_btpv)
x = dPVdy[westerly_idx]
y = -hypvis_btpv[westerly_idx]
nonnan_idx = ~np.isnan(np.log(y))
K0, mu_d_beta2d = curve_fit(lambda x, K0, mu_d_beta2: x*K0/(1+mu_d_beta2*x**2), x, y,
    [y.mean()*x.mean(), 6./50**2])[0]
k, b = stats.linregress(np.log10(x[nonnan_idx]), np.log10(y[nonnan_idx]))[:2]
x1 = np.linspace(np.min(x), np.max(x), 100)

plt.subplot(1,2,1)
# plt.scatter(dPVdy[westerly_idx], -hypvis_btpv[westerly_idx], c='r',label='westerly')
# plt.scatter(dPVdy[~westerly_idx],-hypvis_btpv[~westerly_idx], c='g',label='easterly')
plt.scatter(x, y,label='dPV/dy vs diff flux')
plt.plot(x1, x1*K0/(1+mu_d_beta2d*x1**2), label=r'$K_0$=%.1e, $\mu\beta^{-2}$=%.1e' %(K0, mu_d_beta2d))
plt.plot(x1, x1**(k) * 10**b, label='y=%.1e*x^%1.1f' %(10**b, k))
plt.legend(loc='best',fontsize=9)
plt.subplot(1,2,2)
plt.scatter(x, y, label='dPV/dy vs diff flux')
# plt.scatter(dPVdy[westerly_idx], -hypvis_btpv[westerly_idx], c='r',label='westerly')
# plt.scatter(dPVdy[~westerly_idx],-hypvis_btpv[~westerly_idx], c='g',label='easterly')
plt.plot(x1, x1*K0/(1+mu_d_beta2d*x1**2), label=r'$K_0$=%.1e, $\mu\beta^{-2}$=%.1e' %(K0, mu_d_beta2d))
plt.plot(x1, x1**(k) * 10**b, label='y=%.1e*x^%1.1f' %(10**b, k))
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.legend(loc='best',fontsize=9)
plt.show()

## fit dPV/dy vs. enstrophy forcing
dPVdy = dfdy(lang_btpv)
x = dPVdy[westerly_idx]
y = rm_forc[westerly_idx]
nonnan_idx = ~np.isnan(np.log(y))
F0 = curve_fit(lambda x, F0: F0/x, x, y,
    [y.mean()*x.mean()])[0]
k, b = stats.linregress(np.log10(x[nonnan_idx]), np.log10(y[nonnan_idx]))[:2]
x1 = np.linspace(np.min(x), np.max(x), 100)

plt.subplot(1,2,1)
# plt.scatter(dPVdy[westerly_idx], rm_forc[westerly_idx], c='r',label='westerly')
# plt.scatter(dPVdy[~westerly_idx],rm_forc[~westerly_idx], c='g',label='easterly')
plt.scatter(x, y,label='dPV/dy vs ens forc')
plt.plot(x1, F0/x1, label=r'$F_0$=%.1e' %F0)
plt.plot(x1, x1**(k) * 10**b, label='y=%.1e*x^%1.2f' %(10**b, k))
plt.legend(loc='best',fontsize=9)
plt.subplot(1,2,2)
# plt.scatter(dPVdy[westerly_idx], rm_forc[westerly_idx], c='r',label='westerly')
# plt.scatter(dPVdy[~westerly_idx],rm_forc[~westerly_idx], c='g',label='easterly')
plt.scatter(x, y, label='dPV/dy vs ens forc')
plt.plot(x1, F0/x1, label=r'$F_0$=%.1e' %F0)
plt.plot(x1, x1**(k) * 10**b, label='y=%.1e*x^%1.2f' %(10**b, k))
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.legend(loc='best',fontsize=9)
plt.show()

## save figure
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
import os
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

save_name = 'dPVdy_vs_ensforc_zoom1'
plt.savefig(save_dir + save_name)