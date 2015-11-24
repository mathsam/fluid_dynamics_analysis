import sys
sys.path.append('/home/j1c/py/lib')
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style

filename_prefix = 'Nov5_Sc1.6_drag5e-4'
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
filename = r'%s_seg[0-9]+' %filename_prefix

psif = nc_tools.NetCDFChain(filedir, filename,'psi', last_n_files=2)
psic = qg_transform.get_eddy(qg_transform.real2complex(psif[:]))

bt_psi = np.mean(psic, -1)
tau    = 0.5*(psic[...,0]-psic[...,1])

vor = qg_transform.get_vorticity(bt_psi)
Js_btpsi_vor = qg_transform.grid2spec(qg_transform.jacobian(bt_psi, vor))
k, enstropy_trans = qg_transform.prod_spectrum(vor, Js_btpsi_vor)

Js_btpsi_tau = qg_transform.grid2spec(qg_transform.jacobian(bt_psi, -tau))
k, pe_trans = qg_transform.prod_spectrum(tau, Js_btpsi_tau)

Js_btpsi_bcvor = qg_transform.grid2spec(
    qg_transform.jacobian(bt_psi, qg_transform.get_vorticity(tau)))
k, bc2bts = qg_transform.prod_spectrum(tau, Js_btpsi_bcvor)
#ky, bc2bts_ky = qg_transform.prod_spectrum_meridional(tau, Js_btpsi_bcvor)
##
plt.semilogx(k, -k*enstropy_trans/beta**2, label=r"$-\mathrm{Re}[\zeta_{k}^{'*}J_{k}(\psi',\zeta')]$")
plt.semilogx(k, k*pe_trans,label=r"$-\mathrm{Re}[\tau'_{k}^{*}J_{k}(\psi',\tau')]$")
plt.legend(loc='best')
plt.show()