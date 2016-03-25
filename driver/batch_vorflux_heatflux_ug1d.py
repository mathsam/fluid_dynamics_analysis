import sys
sys.path.append('/home/j1c/py/lib')
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style

filename_prefix = 'Nov4_Sc2.0_drag5e-3'
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix


import os
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

for i_seg in range(20, 41):
    filename = r'%s_seg%d+' %(filename_prefix, i_seg)
    psif = nc_tools.NetCDFChain(filedir, filename,'psi',last_n_files=3)
    psic = qg_transform.real2complex(psif[:])
    
    bt_psi = np.mean(psic, -1)
    u, v = qg_transform.get_velocities(bt_psi)
    
    tau = qg_transform.get_eddy(np.diff(psic, axis=-1)*(-0.5))
    vor = qg_transform.get_eddy(qg_transform.get_vorticity(bt_psi))
    del bt_psi
    vg = np.squeeze(qg_transform.spec2grid(v))
    ug = np.squeeze(qg_transform.spec2grid(u))
    del v, u
    
    taug = np.squeeze(qg_transform.spec2grid(tau))
    vorg = np.squeeze(qg_transform.spec2grid(vor))
    del vor, tau
    heatflux = np.mean(vg*taug, -1)
    vorflux = np.mean(vorg*vg, -1)
    ug1d     = np.mean(ug, -1)
    
    matlab_style.save(save_dir + 'heatflux_vorflux_ug1d_seg%d' %i_seg, 
        heatflux=heatflux, vorflux=vorflux, ug1d=ug1d)