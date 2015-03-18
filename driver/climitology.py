import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np

filename_prefix = 'Jan17_drag_1e-4'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
psi = nc_tools.ncread('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                      '%s_seg(180)' %filename_prefix,'psi')[:]
psi = qg_transform.real2complex(psi)
                      
## zonal mean zonal wind                      

uk, vk = qg_transform.get_velocities(psi)
ug  = qg_transform.spec2grid(uk)

mean_ug = np.mean(np.mean(np.mean(ug, -1), -1), 0)

## zonal mean PV
F    = 1583.14
beta = 3957.9
pvg = qg_transform.spec2grid(qg_transform.get_PV(psi, F))
pvg += qg_transform.get_betay(pvg, beta)

mean_pvg = np.mean(np.mean(np.mean(pvg, -1), -1), 0)