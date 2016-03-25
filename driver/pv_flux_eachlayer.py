import sys
sys.path.append('/home/j1c/py/lib')
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style

filename_prefix = 'Nov4_Sc2.0_drag5e-4'
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
filename = r'%s_seg[0-9]+' %filename_prefix
F = 395.785873603

psif = nc_tools.NetCDFChain(filedir, filename,'psi',last_n_files=1)
psic = qg_transform.real2complex(psif[:])

v = qg_transform.get_velocities(psic)[1]
vg = qg_transform.spec2grid(v)
pvc = qg_transform.get_PV(psic, F)
pvg = qg_transform.spec2grid(pvc)
pvflux = np.mean(np.mean(vg*pvg, 0), 1)

##
import os
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)
    
save_name = 'pv_flux_eachlayer.png'

plt.plot(pvflux)
plt.savefig(save_dir + save_name)
plt.show()