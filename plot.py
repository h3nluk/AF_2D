import numpy as np
import matplotlib.pyplot as plt

dat = np.genfromtxt('f.dat')
plt.imshow(dat)
plt.colorbar()
plt.show()

