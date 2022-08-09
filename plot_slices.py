
import matplotlib.pyplot as plt
import numpy as np


f = np.genfromtxt('f.dat')
f0 = np.genfromtxt('f0.dat')

sizex, sizev = np.shape(f)

fig, ax = plt.subplots(nrows = 2, ncols = 1)

ax[0].plot(f[sizex//2, :], label = 'f')
ax[0].plot(f0[sizex//2, :], label = 'f0')
ax[0].set_xlabel('v')
ax[0].grid()
ax[0].legend(loc = 'best')

ax[1].plot(f[:,sizev//2], label = 'f')
ax[1].plot(f0[:,sizev//2], label = 'f0')
ax[1].set_xlabel('x')
ax[1].grid()
ax[1].legend(loc = 'best')

plt.show()
