import qg_transform
import nc_tools
#import matplotlib.pyplot as plt
import numpy as np
import diag

# Sc = 1.6, F = 1583.14, beta = 3957.86
# Sc = 2.0, F = 1583.14, beta = 3166.29
# Sc = 2.5, F = 1583.14, beta = 2533.03
F    = 1583.14
beta = 2533.03
last_n_files = 2

filename_prefix = 'Jan18_c2.5_drag_1e-2'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
psi = nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                      '%s_seg' %filename_prefix,'psi', last_n_files=last_n_files)[:]
psi = qg_transform.real2complex(psi)

mean_u     = diag.zonal_mean_zonal_wind(psi)
mean_pv    = diag.zonal_mean_PV(psi, F, beta)
mean_dpvdy = diag.zonal_mean_dPVdy(psi, F, beta)

vg  = qg_transform.spec2grid(qg_transform.get_velocities(psi)[1])
pvg = qg_transform.spec2grid(qg_transform.get_PV(psi, F))
bt_vg = np.mean(vg, -1)[...,np.newaxis]
pvflux = np.mean(np.mean(vg*pvg, 0), -2)
btpvflux = np.mean(np.mean(bt_vg*pvg, 0), -2)
mean_EKE = np.mean(np.mean(vg*vg, 0), -2)

## save as npz
import matlab_style as mlab
import os

if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

mlab.save(save_dir + 'zonal_mean_field', mean_u=mean_u, mean_pv=mean_pv, 
                                         mean_dpvdy=mean_dpvdy, pvflux=pvflux,
                                         mean_EKE=mean_EKE, btpvflux=btpvflux)

## plot the figures
import matplotlib
import matplotlib.pyplot as plt
font = {'size'   : 12}
matplotlib.rc('font', **font)

def plot_zonal_field(mean_field, title=None, figsize=(4, 6)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    y = np.linspace(-np.pi, np.pi, mean_field.shape[0])
    ax.plot(mean_field[:,0], y, label='upper')
    ax.plot(mean_field[:,1], y, label='lower')
    ax.set_ylabel('y')
    if title:
        ax.set_title(title)
    return fig
    
fig_u  = plot_zonal_field(mean_u, 'u')
fig_pv = plot_zonal_field(mean_pv, 'pv')
fig_dpvdy = plot_zonal_field(mean_dpvdy, 'dPVdy')
fig_pvflux = plot_zonal_field(pvflux, 'PVflux')
fig_btpvflux = plot_zonal_field(btpvflux, 'barotropic PV flux')
fig_EKE = plot_zonal_field(mean_EKE, 'EKE')

fig_u.savefig(save_dir + 'mean_u.png')
fig_pv.savefig(save_dir + 'mean_PV.png')
fig_dpvdy.savefig(save_dir + 'mean_dPVdy.png')
fig_pvflux.savefig(save_dir + 'PVflux.png')
fig_EKE.savefig(save_dir + 'EKE.png')
fig_btpvflux.savefig( save_dir + 'BT_PVflux.png')