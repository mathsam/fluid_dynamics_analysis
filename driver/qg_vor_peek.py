import sys
sys.path.append('/home/j1c/py/lib')
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style

which_layer = 0
filename_prefix = 'Dec12_kf64_drag1e-4'
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
filename = r'%s_seg239' %filename_prefix
psif = nc_tools.ncread(filedir, filename,'psi')
psi = psif[-1,:,:,:]
psic = qg_transform.real2complex(psi)
vors = qg_transform.get_vorticity(psic)
vorg = qg_transform.spec2grid(vors)
vormap = plt.imshow(vorg[...,which_layer])
plt.colorbar()
vormap.set_cmap('gray')
vormap_std = np.std(vorg[...,which_layer].flatten())
vormap.set_clim(-2*vormap_std, 2*vormap_std)
fig = plt.gcf()
plt.show()

## tracer_y

tracerc = qg_transform.real2complex(nc_tools.ncread(filedir, filename,'tracer_y')[-1])
tracerg = qg_transform.spec2grid(tracerc)
plt.imshow(tracerg[...,which_layer], cmap='gray')
plt.show()

## save figure
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
import os
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

save_name = 'vor_layer=%d_t=%d.png' %(which_layer, psif.total_time_steps)
matlab_style.fig_save(fig, save_dir + save_name)