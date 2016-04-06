import sys
sys.path.append('/home/j1c/py/lib')
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style
# filename_set = ['Dec12_kf16_drag256e-4',
#     'Dec12_kf32_drag256e-4',
#     'Dec12_kf64_drag256e-4']
# filename_set = ['Dec12_kf16_drag8e-4',
#     'Dec12_kf32_drag8e-4',
#     'Dec12_kf64_drag8e-4',
#     'Dec12_kf128_drag8e-4']
exp_drag = 1
exp_year = 2015
# filename_set = ['Mar14_kf16_drag%de-4' %exp_drag,
#     'Mar14_kf32_drag%de-4' %exp_drag,
#     'Mar14_kf64_drag%de-4' %exp_drag,
#     'Mar14_kf128_drag%de-4' %exp_drag]
filename_set = ['Dec12_kf16_drag%de-4' %exp_drag,
    'Dec12_kf32_drag%de-4' %exp_drag,
    'Dec12_kf64_drag%de-4' %exp_drag,
    'Dec12_kf128_drag%de-4' %exp_drag]
kfs = np.array([16, 32, 64, 128])
zonal_energy_ratio = np.zeros(kfs.shape, dtype=float)
num_kf = kfs.size
which_layer = 0


for i, filename_prefix in enumerate(filename_set):
    filedir  = '/archive/Junyi.Chai/QG_exp/%d/%s' %(exp_year, filename_prefix)
    filename = r'%s_seg[0-9]+' %filename_prefix
    psic = qg_transform.real2complex(
        nc_tools.NetCDFChain(
        '/archive/Junyi.Chai/QG_exp/%d/%s' %(exp_year, filename_prefix),
        '%s_seg' %filename_prefix,'psi', last_n_files=1))
    u, v = qg_transform.get_velocities(psic)
    u2 = qg_transform.prod_domain_ave_int(u, u).mean()
    v2 = qg_transform.prod_domain_ave_int(v, v).mean()
    zonal_energy_ratio[i] = u2/(u2+v2)
