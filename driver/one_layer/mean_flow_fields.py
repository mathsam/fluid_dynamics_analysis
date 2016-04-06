import numpy as np
import matlab_style as ml
import sys
sys.path.append('/home/j1c/py/lib')
sys.path.append('/home/j1c/py/lib/spectral_analysis')
import qg_transform
import nc_tools
import diag
import barotropic_spec_diag as bt_sp

filename_prefix = 'Dec12_kf64_drag1e-4'
# kf16 1780, 1796
# kf32 1440, 1470
seg_list = range(1360, 1377) 
filedir  = '/archive/Junyi.Chai/QG_exp/2015/%s' %filename_prefix
save_dir = '/home/j1c/analysis/2015/qg_model/%s/spec_diag/' %filename_prefix

import os
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

mean_psic = diag.func_on_multi_files(filedir,
    filename_prefix, seg_list, np.mean,
    'psi', 0)
ml.save(save_dir + 'mean_psic', mean_psic=mean_psic[0])