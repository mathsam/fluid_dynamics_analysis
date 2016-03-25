import sys
sys.path.append('/home/j1c/py/lib/spectral_analysis')
import qg_transform
import nc_tools
import numpy as np
import matlab_style as ml
import spectral_energy_budget as sp

filename_prefix = 'Nov5_Sc1.6_drag5e-4'
seg_range = range(90, 106) 
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
save_dir = '/home/j1c/analysis/2015/qg_model/%s/btEKEspec/' %filename_prefix

import os
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

for seg in seg_range:
    filename = r'%s_seg%d' %(filename_prefix, seg)
    print "working on %s" %filename
    
    psif = nc_tools.ncread(filedir, filename, 'psi')
    psic = qg_transform.real2complex(psif[:])
    
    BtEKEAnalysis = sp.BarotropicEKE(psic)
    k, eddy_eddy_trans = BtEKEAnalysis.eddy_eddy_transfer()
    k, eddy_mean_trans = BtEKEAnalysis.eddy_mean_transfer()
    k, eddy_eddy_gens  = BtEKEAnalysis.eddy_eddy_gen()
    k, eddy_mean_gens  = BtEKEAnalysis.eddy_mean_gen()

    ml.save(save_dir + 'bt_EKEspec' + str(seg), k=k, eddy_eddy_trans=eddy_eddy_trans,
        eddy_mean_trans=eddy_mean_trans, eddy_eddy_gens=eddy_eddy_gens, eddy_mean_gens=eddy_mean_gens)


##
filename_prefix = 'Nov4_drag5e-4'
seg_range = range(90, 121) 
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
save_dir = '/home/j1c/analysis/2015/qg_model/%s/btEKEspec/' %filename_prefix

eddy_eddy_trans_mean = np.zeros(511)
eddy_mean_trans_mean = np.zeros(511)
eddy_eddy_gens_mean  = np.zeros(511)
eddy_mean_gens_mean  = np.zeros(511)

for seg in seg_range:
    ml.load(save_dir + 'bt_EKEspec' + str(seg) + '.npz')
    eddy_eddy_trans_mean += eddy_eddy_trans
    eddy_mean_trans_mean += eddy_mean_trans
    eddy_eddy_gens_mean  += eddy_eddy_gens
    eddy_mean_gens_mean  += eddy_mean_gens

eddy_eddy_trans_mean  /= len(seg_range)
eddy_mean_trans_mean  /= len(seg_range)
eddy_eddy_gens_mean   /= len(seg_range)
eddy_mean_gens_mean   /= len(seg_range)
#

import matplotlib.pyplot as plt
plt.semilogx(k, k*eddy_eddy_trans_mean, label='eddy-eddy transfer')
plt.semilogx(k, k*eddy_mean_trans_mean, '--b', label='eddy-mean')
plt.semilogx(k, k*eddy_eddy_gens_mean, 'r', label='eddy-eddy gen')
plt.semilogx(k, k*eddy_mean_gens_mean, '--r', label='eddy-mean gen')
plt.xlabel('Wavenumber')
plt.ylabel('Barotropic EKE')
plt.legend(loc='best')
plt.savefig(save_dir + 'btEKEspecBudget.png')
plt.show()