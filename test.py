import PIL
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot([-5, -4, -3, -2 ,-1, 0, -1, -2, -3, -4, -5], [0, 0, 0, 0 , 0, 0, 0, 0, 0, 0, 0], [10, 8, 6, 4, 2, 0, -2, -4 ,-6, -8, -10], label='Satellite-orbit', color='r')
ax.legend()

plt.show()