import numpy as np
import matlab_style as ml
import sys
sys.path.append('/home/j1c/py/lib')
sys.path.append('/home/j1c/py/lib/noboru')
import qg_transform
import nc_tools
import diag
from finite_amp_diag import time_mean_wave_activity_op

filename_prefix = 'Dec12_kf16_drag1e-4'
# kf16 1780, 1801
# kf32 1450, 1471
# kf64 1360, 1381
seg_list = range(1780, 1801) 
filedir  = '/archive/Junyi.Chai/QG_exp/2015/%s' %filename_prefix
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix

import os
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

def PV_lang_mean(psic, beta=16*np.pi):
    vorc = qg_transform.get_vorticity(psic)
    vorg = qg_transform.spec2grid(vorc).squeeze()
    lang_PV, A_PV = time_mean_wave_activity_op(vorg, 1024, zoom=1, beta=beta)
    return lang_PV, A_PV

lang_PV, A_PV = diag.func_on_multi_files(filedir,
    filename_prefix, seg_list, PV_lang_mean, 'psi')

ml.save(save_dir + 'lang_mean_PV', lang_PV=lang_PV, A_PV=A_PV)
