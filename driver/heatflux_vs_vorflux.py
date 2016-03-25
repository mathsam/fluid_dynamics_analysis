import sys
sys.path.append('/home/j1c/py/lib')
import qg_transform
import nc_tools
import matplotlib.pyplot as plt
import numpy as np
import matlab_style

filename_prefix = 'Nov4_Sc2.0_drag5e-4'
filedir  = '/archive/Junyi.Chai/QG_exp/%s' %filename_prefix
filename = r'%s_seg[0-9]+' %filename_prefix


psif = nc_tools.NetCDFChain(filedir, filename,'psi',last_n_files=1)
psic = qg_transform.real2complex(psif[:])

bt_psi = np.mean(psic, -1)
u, v = qg_transform.get_velocities(bt_psi)

tau = np.squeeze(qg_transform.get_eddy(np.diff(psic, axis=-1)*(-0.5)))
vor = qg_transform.get_eddy(qg_transform.get_vorticity(bt_psi))
del bt_psi
vg = np.squeeze(qg_transform.spec2grid(v))
del v
u_ave = np.mean(u, 0)[np.newaxis,...]
ug_ave = np.squeeze(qg_transform.spec2grid(u_ave))
taug = np.squeeze(qg_transform.spec2grid(tau))
vorg = np.squeeze(qg_transform.spec2grid(vor))
del vor, tau

ug1d = np.mean(ug_ave, -1)
lat = np.arange(0, 1024)
#easterly_idx = np.logical_and(np.logical_and(ug1d<0, lat>100), lat<900)
easterly_idx = ug1d<0

heatflux = np.mean(np.mean(vg*taug, 0), -1)
vor_flux = np.mean(np.mean(vorg*vg, 0), -1)
## plot heat flux and vor flux
beta = 791.571747206000 #3957.85873603000
import os
save_dir = '/home/j1c/analysis/2015/qg_model/%s/' %filename_prefix
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)

plt.plot(vor_flux/beta,label='bt_vor flux/beta')
plt.plot(heatflux,label='heat flux')
plt.plot(ug1d/ug1d.std()*heatflux.std(),label='rescaled U')
plt.legend(loc='best')
#plt.savefig(save_dir + 'heatflux_vorflux.png')
plt.show()
##
vorg1d = np.mean(np.mean(vorg, 0), -1)
vorg1d.shape = (1, 1024, 1)
taug1d = np.mean(np.mean(taug, 0), -1)
taug1d.shape = (1, 1024, 1)

##
taug -= taug1d
vorg -= vorg1d

vor_flux = np.mean(np.mean(vorg*vg, 0), -1)
##
v = qg_transform.filter(v, 30, None)
tau = qg_transform.filter(tau, 30, None)
vg = np.squeeze(qg_transform.spec2grid(v))
taug = np.squeeze(qg_transform.spec2grid(tau))
heatflux = np.mean(np.mean(vg*taug, 0), -1)

##
plt.plot(heatflux)
plt.plot(ug1d/ug1d.std()*heatflux.std())
plt.show()
##
plt.plot(vor_flux)
plt.plot(ug1d/ug1d.std()*vor_flux.std())
plt.show()
##
k, heatflux_spec = map(np.squeeze, qg_transform.prod_spectrum(v, tau))
k, v2_spec = map(np.squeeze, qg_transform.prod_spectrum(v, v))

##
#beta = 3957.85873603000 #3957.85873603000# 989.464684007
vorflux3d  = vg[:,easterly_idx,:]*vorg[:,easterly_idx,:]/beta
heatflux3d = vg[:,easterly_idx,:]*taug[:,easterly_idx,:]

print "ration of vorflux/beta/heatflux in easterjet :"
vf_over_hf = -vorflux3d.mean()/heatflux3d.mean()
print vf_over_hf
CUT_OFF = 0.4
vorflux_freq, vorflux_bins = np.histogram(-vorflux3d.flatten()[np.abs(vorflux3d.flatten())<CUT_OFF], 
    4000, [-CUT_OFF, CUT_OFF])
tauflux_freq, tauflux_bins = np.histogram(heatflux3d.flatten()[np.abs(heatflux3d.flatten())<CUT_OFF], 
    4000, [-CUT_OFF, CUT_OFF])

x_vorflux = 0.5*(vorflux_bins[:-1]+vorflux_bins[1:])
x_tauflux = 0.5*(tauflux_bins[:-1]+tauflux_bins[1:])
##
plt.semilogy(x_vorflux, vorflux_freq,'b')
#plt.semilogy(vorflux_bins[1:], vorflux_freq,'--g')
#plt.semilogy(x_mean, y_freq)
plt.semilogy(x_tauflux, tauflux_freq,'r')
plt.xlabel(r"$v'\tau'$, $-v'\zeta'/\beta$")
plt.ylabel('freq')
plt.savefig(save_dir + 'heatflux_vorflux_vs_freq.png')
plt.show()

##
vorflux_netfreq = (vorflux_freq[x_vorflux>0] - vorflux_freq[x_vorflux<0][::-1])
tauflux_netfreq = tauflux_freq[x_tauflux>0] - tauflux_freq[x_tauflux<0][::-1]
vorflux_absfreq = (vorflux_freq[x_vorflux>0] + vorflux_freq[x_vorflux<0][::-1])
tauflux_absfreq = tauflux_freq[x_tauflux>0]  + tauflux_freq[x_tauflux<0][::-1]
x = x_vorflux[x_vorflux>0].copy()
plt.subplot(2,1,1)
plt.plot(x, vorflux_netfreq,'b')
plt.plot(x, tauflux_netfreq,'r')
plt.ylabel('net freq')
plt.subplot(2,1,2)
plt.plot(x, np.cumsum(x*vorflux_netfreq)/np.cumsum(vorflux_absfreq),'b',label='mean vor flux/beta')
plt.plot(x, np.cumsum(x*tauflux_netfreq)/np.cumsum(tauflux_absfreq),'r',label='mean tau flux')
plt.legend(loc='best')
plt.xlabel(r"$v'\tau'$ bin")
plt.savefig(save_dir + 'net_heatflux_vorflux_vs_freq.png')
plt.show()
##
def plot_even(x, y, option=''):
    plt.plot(x[x>0], y[x>0], option)
    plt.plot(-x[x<0], y[x<0], option+'--')

plot_even(x_vorflux, vorflux_freq, 'b')
plot_even(x_tauflux, tauflux_freq, 'r')
plt.show()
##
plt.plot(-vorflux_bins[1:]/beta*0.6*2, vorflux_freq,'b')
plt.plot(-2*5.5e-4 + vorflux_bins[1:]/beta*0.6*2, vorflux_freq,'--b')
plt.plot(x_mean, y_freq)
plt.plot(tauflux_bins[1:]*2, tauflux_freq,'r')
plt.show()
##
some_flux_bins = (vorflux_bins[:-1]+vorflux_bins[1:])*0.5*(-1)
some_flux_freq = vorflux_freq
y_freq = 10**np.arange(7, 0, -0.5)
x_mean = np.zeros_like(y_freq)

for i, y in enumerate(y_freq):
    idx = some_flux_freq >= y
    x_mean[i] = np.sum(some_flux_bins[idx]*some_flux_freq[idx])/np.sum(some_flux_freq[idx])