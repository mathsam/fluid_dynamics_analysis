import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style as mlab
import os

filename_prefix = 'Jan17_drag_1e-4'
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

for i in range(1, 181):
    psi = nc_tools.NetCDFChain('/archive/Junyi.Chai/QG_exp/%s' %filename_prefix,
                        '%s_seg%d.nc' %(filename_prefix, i),'psi', last_n_files=1)
    psi = qg_transform.real2complex(psi)
    k,Ek,EKEk = qg_transform.barotropic_Ek(psi)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog(k, Ek, label='KE')
    ax.loglog(k, EKEk, ls='--', label='EKE')
    most_energetic_k = np.where(EKEk == EKEk.max())[0][0]
    x   = np.arange(10,101)
    p5  = Ek[most_energetic_k]*x**(-5.0)/(x[most_energetic_k]**(-5.0))/5.0
    p53 = Ek[most_energetic_k]*x**(-5.0/3.0)/(x[most_energetic_k]**(-5.0/3.0))/100.0
    ax.loglog(x,p5,  label='$k^{-5}$')
    ax.loglog(x,p53, label='$k^{-5/3}$')
    ax.set_xlabel('Wavenumber')
    ax.set_ylabel('Energy spectrum')
    ax.legend(loc='best')
    #plt.show()

    fig.savefig(save_dir + 'barotropic_Ek_seg%d.png' %i)
    mlab.save(save_dir + 'Ek_seg%d' %i, k=k, Ek=Ek, EKEk=EKEk)