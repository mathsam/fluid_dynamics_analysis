import nc_tools
import qg_transform
import numpy as np
from math import pi

exp_list  = ['Jan17_drag_', 'Jan18_c2_drag_', 'Jan18_c2.5_drag_']
cri_list  = ['1.6', '2.0', '2.5']
drag_list = ['1e0', '1e-1', '1e-2', '1e-3', '5e-4', '1e-4']
arch_dir  = '/archive/Junyi.Chai/QG_exp/'
average_period = 200
U_BAR    = 1.
ROSSBY_R = 2*pi/500.

num_rows = len(drag_list)
num_colmns = len(cri_list)

heat_flux_int = np.zeros((num_rows, num_colmns))

for exp_i, exp_name in enumerate(exp_list):
    for drag_i, drag_name in enumerate(drag_list):
        filedir = arch_dir + exp_name + drag_name
        filename = r'%s_seg[0-9]+' %(exp_name + drag_name)
        psi = nc_tools.NetCDFChain(filedir,
            filename, 'psi', last_n_files=1)[:]
        psi = qg_transform.real2complex(psi)
        tau = 0.5*(psi[...,0] - psi[...,1])[..., np.newaxis]
        bt_psi = np.mean(psi, axis=-1)[..., np.newaxis]
        bt_v = qg_transform.get_velocities(bt_psi)[1]
        heat_flux_int[drag_i, exp_i] = np.mean(qg_transform.prod_domain_ave_int(tau, bt_v))