import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style as mlab
import os

filename_prefix = 'Dec12_kf16_drag8e-4'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
psi = nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/2015/%s' %filename_prefix,
                     '%s_seg' %filename_prefix,'psi', last_n_files=20)
# psi = qg_transform.real2complex(nc_tools.ncread('/archive/Junyi.Chai/QG_exp/2015/%s' %filename_prefix,
#                       '%s_seg40' %filename_prefix, 'psi'))
k,Ek,EKEk = qg_transform.barotropic_Ek(psi)
#

if not os.path.isdir(save_dir):
    os.mkdir(save_dir)
mlab.save(save_dir + 'Ek', k=k, Ek=Ek, EKEk=EKEk)
## plot the barotropic spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(k, Ek, label='KE')
ax.loglog(k, EKEk, ls='--', label='EKE')
most_energetic_k = np.where(EKEk == EKEk.max())[0][0]
x   = np.arange(10,201)
p5  = Ek[most_energetic_k]*x**(-5.0)/(x[most_energetic_k]**(-5.0))/5.0
p53 = Ek[most_energetic_k]*x**(-5.0/3.0)/(x[most_energetic_k]**(-5.0/3.0))/100.0
p3 = Ek[most_energetic_k]*x**(-3.0)/(x[most_energetic_k]**(-3.0))
#p4  = Ek[most_energetic_k]*x**(-4.0)/(x[most_energetic_k]**(-4.0))/100.0
ax.loglog(x,p5,  label='$k^{-5}$')
ax.loglog(x,p53, label='$k^{-5/3}$')
ax.loglog(x,p3,  label='$k^{-3}$')
#ax.loglog(x,p4, label='$k^{-4}$')
ax.set_xlabel('Wavenumber')
ax.set_ylabel('Energy spectrum')
ax.legend(loc='best')
plt.show()
# save figure and spectrum data


fig.savefig(save_dir + 'barotropic_Ek.png')

##
filename_prefix = 'Dec12_kf16_drag64e-4'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
#psi = nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/2015/%s' %filename_prefix,
#                      '%s_seg' %filename_prefix,'psi', last_n_files=1)
psi = qg_transform.real2complex(nc_tools.ncread('/archive/Junyi.Chai/QG_exp/2015/%s' %filename_prefix,
                      '%s_seg25' %filename_prefix, 'psi'))
u, v = qg_transform.get_velocities(psi)
k, u2k = qg_transform.prod_spectrum(u,u)
k, v2k = qg_transform.prod_spectrum(v,v)
plt.loglog(k, u2k, label=r'$u^2$')
plt.loglog(k, v2k, label=r'$v^2$')
plt.legend(loc='best')
fig = plt.gcf()
plt.show()
mlab.fig_save(save_dir + 'u2_v2_spec.png')                      