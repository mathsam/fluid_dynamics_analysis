## calculate eddy scale and jet scale
# ------------------------------------------------------------------------------
# eddy scale is defined as (inverse) centroid of EKE
# jet scale is defined as (inverse) centroid of ZKE
# ------------------------------------------------------------------------------
import os
import numpy as np
import diag
import matplotlib.pyplot as plt
import matlab_style as ml
analysis_dir = '/home/j1c/analysis/2015/qg_model'
save_dir = '/home/j1c/analysis/2015/qg_model/scales'

exp_list  = ['Jan17_drag_', 'Jan18_c2_drag_', 'Jan18_c2.5_drag_']
cri_list  = ['1.6', '2.0', '2.5']
drag_list = ['1e0', '1e-1', '1e-2', '1e-3', '5e-4', '1e-4']

cent_EKE = np.zeros((len(cri_list), len(drag_list)))
cent_ZKE = np.zeros((len(cri_list), len(drag_list)))
ZKE = np.zeros((len(cri_list), len(drag_list)))
EKE = np.zeros((len(cri_list), len(drag_list)))

for i, exp_i in enumerate(exp_list):
    for j, drag_i in enumerate(drag_list):
        exp_name = exp_i + drag_i
        file_Ek = os.path.join(analysis_dir, exp_name, 'Ek.npz')
        Ekzip = np.load(file_Ek)
        k    = Ekzip['k']
        EKEk = Ekzip['EKEk']
        Ek   = Ekzip['Ek']
        ZKEk = Ek - EKEk
        
        cent_EKE[i,j] = diag.centroid(k, EKEk)
        cent_ZKE[i,j] = diag.centroid(k, ZKEk)
        ZKE[i,j] = np.sum(ZKEk)
        EKE[i,j] = np.sum(EKEk)
        
ml.save(os.path.join(save_dir, 'energy_scale.npz'), cent_EKE=cent_EKE,
                                                    cent_ZKE=cent_ZKE,
                                                    ZKE=ZKE,
                                                    EKE=EKE)
## plot the ratio of zonal jet scale over eddy scale
fig = plt.figure()
ax  = fig.add_subplot(111)
ax.plot(drag_list, (cent_EKE/cent_ZKE).T, '--o')
ax.set_xscale('log')
ax.set_xlabel('drag')
ax.set_ylabel('jet scale/eddy scale')
plt.legend(cri_list)
plt.show()