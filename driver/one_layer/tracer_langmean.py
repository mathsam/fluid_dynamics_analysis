import sys
sys.path.append('/home/j1c/py/lib')
sys.path.append('/home/j1c/py/lib/noboru')
from finite_amp_diag import OneLayerWaveActivity, vorflux_1d, dfdy
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style
import os

filename_prefix = 'Dec12_kf64_drag256e-4'
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
filename = r'%s_seg[0-9]+' %filename_prefix

tracer_y = qg_transform.real2complex(
    nc_tools.NetCDFChain(filedir, filename,'tracer_y',last_n_files=1))[...,0]
tracer_yg = qg_transform.spec2grid(tracer_y)
lang_ty, A_ty = time_mean_wave_activity_op(tracer_yg, 1024)

save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

matlab_style.save(save_dir + 'lang_ty', 
    lang_ty = lang_ty,
    A_ty=A_ty)