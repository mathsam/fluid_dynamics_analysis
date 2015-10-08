import qg_transform
import nc_tools
import numpy as np
import diag

# Sc = 1.6, F = 1583.14, beta = 3957.86
# Sc = 2.0, F = 1583.14, beta = 3166.29
# Sc = 2.5, F = 1583.14, beta = 2533.03
F    = 1583.14
beta = 3957.86

filename_prefix = 'Jan17_drag_1e-4'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
psi = nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                      '%s_seg' %filename_prefix,'psi', last_n_files=3)[:]
psi = qg_transform.real2complex(psi)
v = qg_transform.get_velocities(psi)[1]
vg = qg_transform.spec2grid(v)
bt_vg = np.mean(vg, -1)
mean_v2 = np.mean(np.mean(vg**2, -2), 0)

mean_u     = diag.zonal_mean_zonal_wind(psi)
mean_pv    = diag.zonal_mean_PV(psi, F, beta)
mean_dpvdy = diag.zonal_mean_dPVdy(psi, F, beta)

## save as npz
import matlab_style as mlab
import os

if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

mlab.save(save_dir + 'zonal_mean_field', mean_u=mean_u, mean_pv=mean_pv, 
                                         mean_dpvdy=mean_dpvdy)

## plot the figures
import matplotlib.pyplot as plt
load_from_npz = False

if load_from_npz:
    import matlab_style as mlab
    import numpy as np
    filename_prefix = 'Jan17_drag_1e0'
    save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
    mlab.load(save_dir + 'zonal_mean_field.npz')

import matplotlib
font = {'size'   : 12}
matplotlib.rc('font', **font)

def plot_zonal_field(mean_field, title=None, figsize=(4, 6)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    y = np.linspace(-np.pi, np.pi, mean_field.shape[0])
    if len(mean_field.shape) == 2 and mean_field.shape[1] != 1:
        ax.plot(mean_field[:,0], y, label='upper')
        ax.plot(mean_field[:,1], y, label='lower')
    else:
        ax.plot(mean_field, y)
    ax.set_ylabel('y')
    ax.legend(loc='best')
    if title:
        ax.set_title(title)
    return fig
    
fig_u  = plot_zonal_field(mean_u, 'u')
fig_pv = plot_zonal_field(mean_pv, 'pv')
fig_dpvdy   = plot_zonal_field(mean_dpvdy, 'dPVdy')
fig_dpvdybt = plot_zonal_field(np.mean(mean_dpvdy, 1), 'barotropic dPVdy')
fig_v2 = plot_zonal_field(mean_v2, "$v'^2$")

fig_u.savefig(save_dir + 'mean_u.png')
fig_pv.savefig(save_dir + 'mean_PV.png')
fig_dpvdy.savefig(save_dir + 'mean_dPVdy.png')
fig_dpvdybt.savefig(save_dir + 'mean_barotropic_dPVdy.png')
fig_v2.savefig(save_dir + 'v2.png')