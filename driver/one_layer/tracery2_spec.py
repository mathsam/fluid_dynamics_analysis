import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style as mlab
import os

filename_prefix = 'Dec12_kf16_drag8e-4'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
tracer_y = nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                      '%s_seg[0-9]+' %filename_prefix,'tracer_y', last_n_files=10)
k,t2, eddy_t2 = qg_transform.barotropic_spec(tracer_y)

if not os.path.isdir(save_dir):
    os.mkdir(save_dir)
mlab.save(save_dir + 'ty_spec', k=k, t2=t2, eddy_t2=eddy_t2)