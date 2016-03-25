import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style as mlab

# filename_set = ['Dec12_kf16_drag256e-4',
#     'Dec12_kf32_drag256e-4',
#     'Dec12_kf64_drag256e-4',
#     'Dec12_kf128_drag256e-4']
# kfs = np.array([16, 32, 64, 128])
filename_set = ['Dec12_kf16_drag8e-4',
    'Dec12_kf32_drag8e-4',
    'Dec12_kf64_drag8e-4']
kfs = np.array([16, 32, 64])
total_EKE = np.zeros(kfs.size)
total_KE  = np.zeros(kfs.size)
total_v2  = np.zeros(kfs.size)


for i, filename_prefix in enumerate(filename_set):
    save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
    mlab.load(save_dir + 'Ek')
    total_EKE[i] = np.sum(EKEk)
    total_KE[i]  = np.sum(Ek)
    mlab.load(save_dir + 'tracer_spec')
    total_v2[i]  = np.sum(v2_k)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogx(kfs, total_EKE,'-o', label='EKE')
ax.semilogx(kfs, total_KE,'-o', label='KE')
ax.semilogx(kfs, total_v2,'-o', label=r'$v^2$')
ax.set_xlabel('kf')
ax.set_ylabel('Energy')
ax.set_xlim([10, 100])
ax.legend(loc='best')
plt.show()