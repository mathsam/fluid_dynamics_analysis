import sys
sys.path.append('/home/j1c/py/lib')
sys.path.append('/home/j1c/py/lib/noboru')
import qg_transform
import nc_tools
import numpy as np
import matlab_style
from finite_amp_diag import OneLayerWaveActivity

filename_prefix = 'Dec12_kf16_drag8e-4_m4'
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
filename = r'%s_seg153.nc' %filename_prefix
psic_bt  = qg_transform.real2complex(
    nc_tools.ncread(filedir, filename,'psi'))[...,0]
rm_forcing = qg_transform.real2complex(
    nc_tools.ncread(filedir, filename,'rm_forcing'))

params = {'beta': 50.2654824574000, 
          'bot_drag' : 4*8e-4}

fwd = OneLayerWaveActivity(psic_bt, params)
fwd.cal_barotropic(rm_forcing=rm_forcing)

u, v = qg_transform.get_velocities(psic_bt)

u_ave = np.mean(u, 0)[np.newaxis,...]
ug_ave = np.squeeze(qg_transform.spec2grid(u_ave))
ug1d = np.mean(ug_ave, -1)

lat_grid = np.arange(0, 1024)
westerly_idx = ug1d>0

vorflux = vorflux_1d(psic_bt)
#
hypvis_btpv = fwd.hypvis_btpv.copy()
drag_btpv = fwd.drag_btpv.copy()
lang_btpv = fwd.lang_btpv.copy()
rm_forc   = fwd.rm_forc.copy()
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
import os
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)
matlab_style.save(save_dir + 'noboru_budget_seg153', 
    hypvis_btpv = hypvis_btpv, drag_btpv = drag_btpv,
    lang_btpv = lang_btpv,
    ug1d = ug1d,
    rm_forc=rm_forc,
    vorflux=vorflux)