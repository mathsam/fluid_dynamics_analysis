# drag=256e-4
epsilon = 0.0157913670417
beta    = 50.2654824574
kfs = np.array([16, 32, 64, 128])
tracer_flux = np.array([0.001488, 0.003888, 0.006619, 0.008714])
plt.loglog(kfs, tracer_flux, '-o', label='r=256e-4')
plt.loglog(kfs, beta**(-2)*epsilon*kfs**2, label=r'$\beta^{-2}k_f^2\varepsilon$')
plt.legend(loc='best')
#plt.xlim([10, 200])
#plt.xlabel(r'$k_f$')
#plt.ylabel('Tracer diffusivity')

#plt.show()

# drag=8e-4
epsilon = 4.934802200540000E-004
beta    = 50.2654824574
kfs = np.array([16, 32, 64])
tracer_flux = np.array([0.000027, 0.000120, 0.000338])
plt.loglog(kfs, tracer_flux, '-o', label='r=8e-4')
plt.loglog(kfs, beta**(-2)*epsilon*kfs**2, label=r'$\beta^{-2}k_f^2\varepsilon$')
plt.legend(loc='best')
plt.xlim([10, 200])
plt.xlabel(r'$k_f$')
plt.ylabel('Tracer diffusivity')

plt.show()