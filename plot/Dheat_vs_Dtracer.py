plt.scatter(D_heat, tracerflux, label=r'$D_{heat}$ vs $D_{tracer}$',c='r')
x = np.linspace(np.min(D_heat), np.max(D_heat), 10)
y = 2.463*x
plt.plot(x, y, label='y = 2.46x')
plt.legend(loc='best')
plt.xlabel(r'$D_{heat}$')
plt.ylabel(r'$D_{tracer}$')
plt.show()

##
plt.scatter(D_heat, tracerflux, label=r'$D_{heat}$ vs $D_{tracer}$',c='r')
x = np.linspace(np.min(D_heat), np.max(D_heat), 10)
y = 2.9* x**(1.035)
plt.plot(x, y, label='$y = 2.9x^{1.035}$')
plt.xlabel(r'$D_{heat}$')
plt.ylabel(r'$D_{tracer}$')
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.legend(loc='lower right')
plt.show()