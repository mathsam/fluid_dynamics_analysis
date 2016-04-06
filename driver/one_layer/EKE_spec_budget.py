import numpy as np
import matlab_style as ml
import sys
sys.path.append('/home/j1c/py/lib')
sys.path.append('/home/j1c/py/lib/spectral_analysis')
import qg_transform
import nc_tools
import diag
import barotropic_spec_diag as bt_sp

filename_prefix = 'Dec12_kf16_drag64e-4'
seg_list = range(70, 91) 
filedir  = '/archive/Junyi.Chai/QG_exp/2015/%s' %filename_prefix
save_dir = '/home/j1c/analysis/2015/qg_model/%s/spec_diag/' %filename_prefix

import os
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

k, eddy_eddy, eddy_mean, total = diag.func_on_multi_files(filedir,
    filename_prefix, seg_list, bt_sp.eddy_mean_energy_transfer,
    'psi')
ml.save(save_dir + 'spectral_EKE_budget', 
    k=k, eddy_eddy=eddy_eddy, eddy_mean=eddy_mean, total=total)
