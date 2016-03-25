import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style as mlab

# filename_set = ['Dec12_kf16_drag256e-4',
#     'Dec12_kf32_drag256e-4',
#     'Dec12_kf64_drag256e-4']
# label_set = ['kf=16','kf=32','kf=64','kf=128']
filename_set = ['Dec12_kf16_drag8e-4',
    'Dec12_kf32_drag8e-4',
    'Dec12_kf64_drag8e-4']
label_set = ['kf=16','kf=32','kf=64']
fig = plt.figure()
ax = fig.add_subplot(111)
for filename_prefix, label in zip(filename_set, label_set):
    save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
    mlab.load(save_dir + 'ty_spec')
    ax.loglog(k, eddy_t2, label=label)


most_energetic_k = np.where(eddy_t2 == eddy_t2.max())[0][0]
x   = np.arange(5,200)
p53 = eddy_t2[most_energetic_k]*x**(-5.0/3.0)/(x[most_energetic_k]**(-5.0/3.0))
ax.loglog(x,p53,'k', label='$k^{-5/3}$')
ax.set_xlabel('Wavenumber')
ax.set_ylabel(r'Eddy tracer spec')
ax.legend(loc='best')
plt.show()