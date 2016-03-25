import qg_transform
import nc_tools
import numpy as np
import diag

# Sc = 1.6, F = 1583.14, beta = 3957.86
# Sc = 2.0, F = 1583.14, beta = 3166.29
# Sc = 2.5, F = 1583.14, beta = 2533.03
#F    = 1583.14
#beta = 2533.03

filename_prefix = 'Dec12_kf16_drag8e-4_m16'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
psi = nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                      '%s_seg' %filename_prefix,'psi', last_n_files=2)[...,0]
psi = qg_transform.real2complex(psi)
t = nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                      '%s_seg' %filename_prefix,'time', last_n_files=2)

u = qg_transform.get_velocities(psi)[0]
ug = qg_transform.spec2grid(u)
del u, psi
zonalmean_ug = np.mean(ug, -1)