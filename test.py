import PIL
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import math


x = [0, 1, 2, 3, 4, 5]
G = 6.673 * (10 ** (-11))
M = 5.972 * (10 ** 24)
print([-1 * G * M * i for i in x[:3]])
print(round((np.linalg.norm(x[:3])), 4))