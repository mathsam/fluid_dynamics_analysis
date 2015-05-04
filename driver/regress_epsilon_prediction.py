import numpy as np

beta = np.array([3957.86, 3166.29, 2533.03])
beta.shape = (1,3)
beta = np.tile(beta, (5,1))

kf = np.array([1e-1, 1e-2, 1e-3, 5e-4, 1e-4])
kf.shape = (5,1)
kf = np.tile(kf, (1, 3))

F = 1583.14
L = np.pi*2
##
epsilon = gen_rate.flatten()
beta = beta.flatten()
kf = kf.flatten()

X = np.zeros((epsilon.size, 2))
X[:,0] = - epsilon**(9./5.)/beta/L/np.sqrt((np.sqrt(F)*kf))
X[:,1] =  4*F*beta**(-0.8)
Y = epsilon**(0.4)

regress = np.linalg.lstsq(X, Y)

##
import functools
import scipy.optimize
def epsilon_eq(x, ibeta, ikf, A, B):
    L = 2*np.pi
    return x**(0.4) + B*x**(0.9)/ibeta/L/np.sqrt(np.sqrt(F)*ikf) - A*4.*F*ibeta**(-0.8)
    
def solve_epsilon(param):
    A = param[0]
    B = param[1]
    pred_epsilon = np.zeros_like(epsilon)
    for i in range(0, pred_epsilon.size):
        curr_epsilon_eq = functools.partial(epsilon_eq, ibeta=beta[i],
            ikf=kf[i], A=A, B=B)
        pred_epsilon[i] = scipy.optimize.newton_krylov(curr_epsilon_eq, 50.)
    return pred_epsilon

def residual(param):
    return solve_epsilon(param) - epsilon

opt_params = scipy.optimize.leastsq(residual, np.array([10., 1.]))
print opt_params

## plot results predicted epsilon vs. simulation epsilon 
myepsilon = solve_epsilon(opt_params[0])
myepsilon.shape = (5, 3)
plt.loglog(gen_rate[:,0], myepsilon[:,0], '--o', label='Sc=1.6')
plt.loglog(gen_rate[:,1], myepsilon[:,1], '--o', label='Sc=2.0')
plt.loglog(gen_rate[:,2], myepsilon[:,2], '--o', label='Sc=2.5')
plt.legend(loc='best')
plt.xlabel(r'predicted $\epsilon$')
plt.ylabel(r'simulation $\epsilon$')
plt.axis('equal')
plt.show()

## plot results kf vs epsilon
x = np.array([1e-1, 1e-2, 1e-3, 5e-4, 1e-4])
plt.loglog(x, myepsilon[:,0], '--ro', label='Sc=1.6 prediction')
plt.loglog(x, gen_rate[:,0], '-rs', label='Sc=1.6 simulation')
plt.loglog(x, myepsilon[:,1], '--bo', label='Sc=2.0 prediction')
plt.loglog(x, gen_rate[:,1], '-bs', label='Sc=2.0 simulation')
plt.loglog(x, myepsilon[:,2], '--go', label='Sc=2.5 prediction')
plt.loglog(x, gen_rate[:,2], '-gs', label='Sc=2.5 simulation')
plt.xlabel('kf')
plt.ylabel(r'$\epsilon$')
plt.legend(loc='best')
plt.show()