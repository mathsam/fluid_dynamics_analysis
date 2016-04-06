import sys
sys.path.append('/home/j1c/py/lib')
sys.path.append('/home/j1c/py/lib/spectral_analysis')
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style as ml
import barotropic_spec_diag as bt_sp

filename_prefix = 'Mar27QL_kf8_drag64e-4'
seg_range = range(19, 20) 

filedir  = '/archive/Junyi.Chai/QG_exp/2016/%s' %filename_prefix
save_dir = '/home/j1c/analysis/2016/qg_model/%s/spec_diag/' %filename_prefix

import os
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

for seg in seg_range:
    filename = r'%s_seg%d' %(filename_prefix, seg)
    print "working on %s" %filename
    psic = qg_transform.real2complex(nc_tools.ncread(filedir, filename,'psi'))
    tracerc = qg_transform.real2complex(nc_tools.ncread(filedir, filename,'tracer_y'))
    #k, trans_k = bt_sp.nonlinear_transfer(psic, tracerc)
    vorc = qg_transform.get_vorticity(psic)
    k, ens_eddy_eddy_trans, ens_eddy_mean_trans, ens_total_trans = bt_sp.eddy_mean_transfer(
        psic, vorc)
    k, tra_eddy_eddy_trans, tra_eddy_mean_trans, tra_total_trans = bt_sp.eddy_mean_transfer(
        psic, tracerc)
    ml.save(save_dir + 'ens_tracer_flux' + str(seg), k=k, 
        ens_eddy_eddy_trans=ens_eddy_eddy_trans, 
        ens_eddy_mean_trans=ens_eddy_mean_trans, 
        ens_total_trans=ens_total_trans,
        tra_eddy_eddy_trans=tra_eddy_eddy_trans, 
        tra_eddy_mean_trans=tra_eddy_mean_trans, 
        tra_total_trans=tra_total_trans)

##
# kf128: 140~144
# kf64: 245~250
# kf32: 349~350
# kf16: 399~401
#filename_prefix = 'Feb24_kf128_drag256e-4'
#seg_range = range(29, 31) 
#filedir  = '/archive/Junyi.Chai/QG_exp/2016/%s' %filename_prefix
#save_dir = '/home/j1c/analysis/2016/qg_model/%s/spec_diag/' %filename_prefix

trans_k_mean = np.zeros(511)

for seg in seg_range:
    ml.load(save_dir + 'tracer_var_flux' + str(seg) + '.npz')
    trans_k_mean += trans_k.flatten()

trans_k_mean  /= len(seg_range)

##
plt.plot(k, np.cumsum(trans_k_kf16)/np.cumsum(trans_k_kf16)[100:250].mean(), label='kf=16')
plt.plot(k, np.cumsum(trans_k_kf32)/np.cumsum(trans_k_kf32)[100:250].mean(), label='kf=32')
plt.plot(k, np.cumsum(trans_k_kf64)/np.cumsum(trans_k_kf64)[100:250].mean(), label='kf=64')
plt.plot(k, np.cumsum(trans_k_kf128)/np.cumsum(trans_k_kf128)[100:250].mean(), label='kf=128')
plt.xlabel('Wavenumber')
plt.ylabel('Normalized tracer variance flux')
plt.legend(loc='best')
plt.show()