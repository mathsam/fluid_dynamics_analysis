#-------------------------------------------------------------------------------
# plot the heat flux spectrum at 3 different times
import matplotlib.pyplot as plt
import numpy as np

data_dir = '/home/j1c/analysis/2015/qg_model/Jan17_drag_1e-4/heatflux/'

data1 = np.load(data_dir + 'flux_spec_seg1.npz')
data2 = np.load(data_dir + 'flux_spec_seg90.npz')
data3 = np.load(data_dir + 'flux_spec_last10ave.npz')

k = data1['k']

fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(k, data1['tau2'], 'b--', label=r"$\tau'^2$, t1")
ax.loglog(k, data1['tau_flux'], 'b-', label=r"$\tau'v$, t1")
ax.loglog(k, data2['tau2'], 'g--', label=r"$\tau'^2$, t2")
ax.loglog(k, data2['tau_flux'], 'g-', label=r"$\tau'v'$, t2")
ax.loglog(k, data3['tau2'], 'k--', label=r"$\tau'^2$, t3")
ax.loglog(k, data3['tau_flux'], 'k-', label=r"$\tau'v'$, t3")

EKEk  = data3['tau2']
Ek = data3['tau2']
most_energetic_k = np.where(EKEk == EKEk.max())[0][0]
x   = np.arange(10,101)
p53 = Ek[most_energetic_k]*x**(-5.0/3.0)/(x[most_energetic_k]**(-5.0/3.0))
ax.loglog(x,p53, label='$k^{-5/3}$')
ax.set_xlabel('Wavenumber')
ax.set_ylabel('Spectrum')
ax.legend(loc='best')
plt.show()

## compare correlation between tau and v at 3 times
corr1 = data1['tau_flux']/np.sqrt(data1['tau2']*data1['v2'])
corr2 = data2['tau_flux']/np.sqrt(data2['tau2']*data2['v2'])
corr3 = data3['tau_flux']/np.sqrt(data3['tau2']*data3['v2'])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogx(k, corr1, label='t1')
ax.semilogx(k, corr2, label='t2')
ax.semilogx(k, corr3, label='t3')
ax.set_xlabel('Wavenumber')
ax.set_ylabel(r"corr($\tau'$, $v'$)")
ax.legend(loc='best')
plt.show()