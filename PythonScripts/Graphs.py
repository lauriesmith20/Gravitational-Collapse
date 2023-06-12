import matplotlib.pyplot as plt
import numpy as np
import random

ParticleSize = np.array([1,2,3,4,5,6,7])
InitialDensity = np.array([10,20,30,40])
Times = np.array([[1,2,3,4],[2,3,4,5],[2,3,4,5],[10,7,6,5],[1,3,4,5],[8,8,8,8],[1,1,1,1]])


fig = plt.figure()
X,Y = np.meshgrid(InitialDensity, ParticleSize )

ax = fig.add_subplot(111, projection='3d')

Z = Times

ax.plot_surface(X, Y , Times, cmap = 'viridis')
plt.show()

#plt.plot(InitialDensity,Times[6,:])
#plt.yscale("log")
#plt.xscale("log")