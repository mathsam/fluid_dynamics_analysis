import sys
sys.path.append('/home/j1c/py/lib/spectral_analysis')
import qg_transform
import nc_tools
import numpy as np
import matlab_style as ml
import spectral_energy_budget as sp

filename_prefix = 'Nov5_Sc1.6_drag5e-1'
seg_range = range(4, 11)
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
save_dir = '/home/j1c/analysis/2015/qg_model/%s/btKEspec/' %filename_prefix

import os
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

for seg in seg_range:
    filename = r'%s_seg%d' %(filename_prefix, seg)
    print "working on %s" %filename
    
    psif = nc_tools.ncread(filedir, filename, 'psi')
    psic = qg_transform.real2complex(psif[:])
    
    BtKEAnalysis = sp.BarotropicKE(psic)
    k, NLTrans = BtKEAnalysis.nonlinear_transfer()
    k, gens    = BtKEAnalysis.ke_generation()
    k, gens_ms = BtKEAnalysis.ke_generation_shear(U=1)
    k, frics   = BtKEAnalysis.friction_diss(0.039788735773) #19.8943678865
    #k, hypdis  = BtKEAnalysis.hyper_diss()
    ml.save(save_dir + 'bt_KEspec' + str(seg), k=k, NLTrans=NLTrans,
        gens=gens, gens_ms=gens_ms, frics=frics)

##
filename_prefix = 'Nov5_Sc1.6_drag5e-4'
seg_range = range(90, 106)
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
save_dir = '/home/j1c/analysis/2015/qg_model/%s/btKEspec/' %filename_prefix

NLTrans_mean = np.zeros(511)
gens_mean    = np.zeros(511)
gens_ms_mean = np.zeros(511)
frics_mean   = np.zeros(511)

for seg in seg_range:
    ml.load(save_dir + 'bt_KEspec' + str(seg) + '.npz')
    NLTrans_mean += NLTrans
    gens_mean    += gens
    gens_ms_mean += gens_ms
    frics_mean   += frics

NLTrans_mean /= len(seg_range)
gens_mean    /= len(seg_range)
gens_ms_mean /= len(seg_range)
frics_mean   /= len(seg_range)
#

import matplotlib.pyplot as plt
plt.semilogx(k, k*NLTrans_mean, label='Nonlinear transfer')
plt.semilogx(k, k*gens_mean, '--r', label='Gen')
plt.semilogx(k, k*gens_ms_mean, 'r', label='Mean shear gen')
plt.semilogx(k, k*frics_mean, 'k', label='Frictional dissp')
plt.xlabel('Wavenumber')
plt.ylabel('Barotropic KE')
plt.legend(loc='best')
plt.savefig(save_dir + 'btKEspecBudget.png')
plt.show()
