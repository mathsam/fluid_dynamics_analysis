"""Calculate the zonal mean flow and add a perburtation to it to create a restart.nc"""
import sys
sys.path.append('/home/j1c/py/lib')
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np

filename_prefix = 'Dec12_kf16_drag1e-4'
filedir  = '/archive/Junyi.Chai/QG_exp/2015/%s' %filename_prefix
filename = r'%s_seg1275' %filename_prefix
psic = qg_transform.real2complex(nc_tools.ncread(filedir, filename,'psi'))
psic = np.mean(psic, 0)
psic_zonal = psic - qg_transform.filter(psic, None, None, True)
uc = qg_transform.get_velocities(psic_zonal)[0]
ug = qg_transform.spec2grid(uc)
# vorc_zonal = qg_transform.get_vorticity(psic_zonal)
# vorg_zonal = qg_transform.spec2grid(vorc_zonal)
# plt.imshow(vorg_zonal[...,0])
# plt.show()

## perturbation is added to spectral ring band [K_MIN, K_MAX]
AMPLITUDE = 0.5
K_MIN = 16
K_MAX = 18
nc_name = 'restart_amp5e-1.nc'
np.random.seed(1)
perturbation = AMPLITUDE/1512*(np.random.normal(0, 1.0, psic.shape) 
                             + 1j*np.random.normal(0, 1.0, psic.shape))
psic_perturbed = qg_transform.filter(perturbation, 14, 18, True) + psic_zonal

f = nc_tools.CreateRestartNC(nc_name)
psic_perturbed.shape = (1, 512, 1023, 1)
f.write_var('psi',psic_perturbed)
del f