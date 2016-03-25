import sys
sys.path.append('/home/j1c/py/lib')
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style
# filename_set = ['Dec12_kf16_drag256e-4',
#     'Dec12_kf32_drag256e-4',
#     'Dec12_kf64_drag256e-4']
filename_set = ['Dec12_kf16_drag8e-4',
    'Dec12_kf32_drag8e-4',
    'Dec12_kf64_drag8e-4']
kfs = np.array([16, 32, 64])
num_kf = kfs.size
which_layer = 0

for i, filename_prefix in enumerate(filename_set):
    filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
    filename = r'%s_seg[0-9]+' %filename_prefix
    psic = qg_transform.real2complex(
        nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
        '%s_seg' %filename_prefix,'psi', last_n_files=1)[-1])
    vors = qg_transform.get_vorticity(psic)
    vorg = qg_transform.spec2grid(vors)
    tracerc = qg_transform.real2complex(
        nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
        '%s_seg' %filename_prefix,'tracer_y', last_n_files=1)[-1])
    tracerg = qg_transform.spec2grid(tracerc)
    plt.subplot(num_kf, 2, 2*i+1)
    plt.axis('off')
    vormap = plt.imshow(vorg[...,which_layer])
    plt.colorbar()
    vormap.set_cmap('gray')
    vormap_std = vorg[...,which_layer].std()
    vormap.set_clim(-2*vormap_std, 2*vormap_std)
    plt.subplot(num_kf, 2, 2*i+2)
    plt.axis('off')
    vormap = plt.imshow(tracerg[...,which_layer])
    plt.colorbar()
    vormap.set_cmap('gray')
    vormap_std = tracerg[...,which_layer].std()
    vormap.set_clim(-2*vormap_std, 2*vormap_std)

fig = plt.gcf()

plt.show()

## tracer_y

tracerc = qg_transform.real2complex(nc_tools.ncread(filedir, filename,'tracer_y')[-1])
tracerg = qg_transform.spec2grid(tracerc)
plt.imshow(tracerg[...,which_layer], cmap='gray')
plt.show()