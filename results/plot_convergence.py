
import numpy as np
import matplotlib.pyplot as plt

# ~ fname = 'error2D_ax=1ay=0.txt'
# ~ fname = 'error2D_ax=1ay=1.txt'
fname = 'error2D_MaxNorm_ax=1ay=0.txt'
# ~ fname = 'error2D_MaxNorm_ax=1ay=1.txt'

dat = np.genfromtxt(fname)

def order(x,A,N):
	return A * x**N

Nxs = dat[:,0]
Nvs = dat[:,1]
errs = np.abs(1. - dat[:,2])

N = Nxs * Nvs

plt.plot(Nxs, errs, marker = 'v', color = 'red', label = 'AF')
plt.plot(Nxs, order(Nxs, 100., -2), color = 'grey', linestyle = ':', label = 'N^-2')
plt.plot(Nxs, order(Nxs, 1000., -3), color = 'black', linestyle = ':', label = 'N^-3')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Nx, Nv')
plt.ylabel('err')
plt.legend(loc = 'best')
plt.grid()
plt.show()
