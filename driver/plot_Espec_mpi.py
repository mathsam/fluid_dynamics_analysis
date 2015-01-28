import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np

psi = nc_tools.ncread('/archive/Junyi.Chai/QG_exp/Jan17_drag_1e-3','Jan17_drag_1e-3_seg(3[0-9]|40)','psi')
k,Ek,EKEk = qg_transform.barotropic_Ek(psi)

## plot the barotropic spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(k, Ek, label='KE')
ax.loglog(k, EKEk, ls='--', label='EKE')
x   = np.arange(10,101)
p5  = Ek[9]*x**(-5.0)/(x[9]**(-5.0))/5.0
p53 = Ek[9]*x**(-5.0/3.0)/(x[9]**(-5.0/3.0))/100.0
ax.loglog(x,p5,  label='k^{-5}')
ax.loglog(x,p53, label='k^{-5/3}')
ax.set_xlabel('Wavenumber')
ax.set_ylabel('Energy spectrum')
ax.legend()
plt.show()
## save figure and spectrum data
import matlab_style as mlab
save_dir = '/home/j1c/analysis/2015/qg_model/Jan17_drag_1e-3/'
fig.savefig(save_dir + 'barotropic_Ek.png')
mlab.save(save_dir + 'Ek', k=k, Ek=Ek, EKEk=EKEk)
