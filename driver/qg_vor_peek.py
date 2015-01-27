import sys
sys.path.append('/home/j1c/py/lib')
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np

which_layer = 1
psi = nc_tools.ncread('/archive/Junyi.Chai/QG_exp/Jan21Reso2x_c2.5_drag_1e-3','Jan21Reso2x_c2.5_drag_1e-3_seg[0-9]+','psi')
psi = psi[-1,:,:,:]
psic = qg_transform.real2complex(psi)
vors = qg_transform.get_vorticity(psic)
vorg = qg_transform.spec2grid(vors)
vormap = plt.imshow(vorg[...,which_layer])
plt.colorbar()
vormap.set_cmap('gray')
vormap_std = np.std(vorg[...,which_layer].flatten())
vormap.set_clim(-2*vormap_std, 2*vormap_std)
plt.show()