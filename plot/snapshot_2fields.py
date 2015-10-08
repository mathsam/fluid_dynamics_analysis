#--------------- Make a 2 panel plot --------------------
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np

filename_prefix = 'Jan17_drag_1e-4'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
psi = nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                      '%s_seg' %filename_prefix,'psi', last_n_files=1)
psi = qg_transform.real2complex(psi)
psi = psi[-1]
bt_psi = 0.5*(psi[...,0] + psi[...,1])[...,np.newaxis]
tau = 0.5*(psi[...,0] - psi[...,1])[...,np.newaxis]
bt_vor = qg_transform.get_vorticity(bt_psi)

bt_vorg = qg_transform.spec2grid(bt_vor)
taug = qg_transform.spec2grid(tau)

eddy_bt_vor = qg_transform.filter(bt_psi, None, None, True)
eddy_tau    = qg_transform.filter(tau, None, None, True)
eddy_bt_vorg = qg_transform.spec2grid(eddy_bt_vor)
eddy_taug    = qg_transform.spec2grid(eddy_tau)
bt_vorg = qg_transform.spec2grid(bt_vor)

## plot eddy fields
field1 = eddy_bt_vorg[...,0]
field2 = eddy_taug[...,0]

plt.subplot(121)
plt.imshow(field1, cmap='gray', interpolation='none')
c_range = np.std(eddy_bt_vorg[...,0].flatten())
plt.clim([-3.*c_range, 3.*c_range])
plt.title(r"$\zeta'$")

plt.subplot(122)
plt.imshow(field2, cmap='gray', interpolation='none')
c_range = np.std(eddy_taug[...,0].flatten())
plt.clim([-3.*c_range, 3.*c_range])
plt.title(r"$\tau'$")
plt.savefig(save_dir + 'snapshot_eddy_vor_tau.png')
plt.show()

## plot total fields
field1 = bt_vorg[...,0]
field2 = taug[...,0]

plt.subplot(121)
plt.imshow(field1, cmap='gray', interpolation='none')
c_range = np.std(eddy_bt_vorg[...,0].flatten())
#plt.clim([-3.*c_range, 3.*c_range])
plt.title(r"$\zeta$")

plt.subplot(122)
plt.imshow(field2, cmap='gray', interpolation='none')
c_range = np.std(eddy_taug[...,0].flatten())
#plt.clim([-3.*c_range, 3.*c_range])
plt.title(r"$\tau$")

plt.savefig(save_dir + 'snapshot_total_vor_tau.png')
plt.show()

##plot zoom-in fields

# plot eddy fields
field1 = eddy_bt_vorg[0:512,0:512,0]
field2 = eddy_taug[0:512,0:512,0]

plt.subplot(121)
plt.imshow(field1, cmap='gray', interpolation='none')
c_range = np.std(eddy_bt_vorg[...,0].flatten())
plt.clim([-3.*c_range, 3.*c_range])
plt.title(r"$\zeta'$")

plt.subplot(122)
plt.imshow(field2, cmap='gray', interpolation='none')
c_range = np.std(eddy_taug[...,0].flatten())
plt.clim([-3.*c_range, 3.*c_range])
plt.title(r"$\tau'$")
plt.savefig(save_dir + 'snapshot_zoomin_eddy_vor_tau.png')
plt.show()

# plot total fields
field1 = bt_vorg[0:512,0:512,0]
field2 = taug[0:512,0:512,0]

plt.subplot(121)
plt.imshow(field1, cmap='gray', interpolation='none')
c_range = np.std(eddy_bt_vorg[...,0].flatten())
#plt.clim([-3.*c_range, 3.*c_range])
plt.title(r"$\zeta$")

plt.subplot(122)
plt.imshow(field2, cmap='gray', interpolation='none')
c_range = np.std(eddy_taug[...,0].flatten())
#plt.clim([-3.*c_range, 3.*c_range])
plt.title(r"$\tau$")

plt.savefig(save_dir + 'snapshot_zoomin_total_vor_tau.png')
plt.show()