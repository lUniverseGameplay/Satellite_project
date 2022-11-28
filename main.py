import PIL
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from orbit import get_orbit
from scipy.integrate import odeint


class Earth():
  def __init__(self, ax):
    self.G = 6.673 * (10 ** (-11))
    self.M = 5.972 * (10 ** 24)
    self.radius = 6371000
    self.Oz = [0, 0, 1]

    self.stratosphere_radius = self.radius + 50000
    
    self.draw_me(ax)
    self.draw_stratosphere(ax)
  
  def draw_me(self, ax):
    # Open Image in PIL
    bm = PIL.Image.open('data/earthicefreesm.jpg') 
    bm = np.array(bm.resize([d * 2 // 1 for d in bm.size]))/256

    # Create Earth with texture
    lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180 
    lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180

    x = self.radius * np.outer(np.cos(lons), np.cos(lats)).T
    y = self.radius * np.outer(np.sin(lons), np.cos(lats)).T
    z = self.radius * np.outer(np.ones(np.size(lons)), np.sin(lats)).T

    # Plot the Earth
    ax.plot_surface(x, y, z, rstride=4, cstride=4, facecolors = bm)
  

  def draw_stratosphere(self, ax):
    # Create Erath's stratosphere
    u_stratosphere = np.linspace(0, 2 * np.pi, 100)
    v_stratosphere = np.linspace(0, np.pi, 100)

    x_stratosphere = self.stratosphere_radius * np.outer(np.cos(u_stratosphere), np.sin(v_stratosphere))
    y_stratosphere = self.stratosphere_radius * np.outer(np.sin(u_stratosphere), np.sin(v_stratosphere))
    z_stratosphere = self.stratosphere_radius * np.outer(np.ones(np.size(u_stratosphere)), np.cos(v_stratosphere))

    # Plot the stratosphere
    ax.plot_wireframe(x_stratosphere, y_stratosphere, z_stratosphere, linewidth=1, alpha=0.3)

  def return_data(self):
    return [self.G, self.M, self.radius, self.stratosphere_radius, self.Oz]


class Satellite():
  def __init__(self, ax, E_data):
    self.mass = 457 * (10 ** 3)
    self.height = 10 ** 6
    self.velocity = 7654
    self.time = 90 * 60
    # self.North = float(input()) * np.pi / 180
    # self.East = float(input())  * np.pi / 180
    self.North = -43.07 * np.pi / 180
    self.East = -61.5  * np.pi / 180

    self.init_pos = [np.cos(self.North) * np.cos(self.East),
    np.cos(self.North) * np.sin(self.East),
    np.sin(self.North)]

    self.orbit_norm = get_orbit(self.init_pos, False)

    self.tau = np.cross(self.orbit_norm, self.init_pos) + [0.5 * i for i in self.orbit_norm]

    r0 = np.asarray([i * (E_data[2] + self.height) for i in self.init_pos])
    v0 = self.tau * self.velocity

    tspan = np.linspace(0, 3 * self.time, 10 ** 5)
    x0 = [r0, v0]
    
    self.odefun = lambda t, x: [np.asarray(x[1]), np.asarray([-1 * E_data[0] * E_data[1] * i / ((np.linalg.norm(i)) ** 3) for i in x[0]])]

    #print(self.odefun(tspan[0], x0))

    [t, x] = odeint(self.odefun, tspan, x0)

    print([t, x])

    self.draw_self_orbit(ax, E_data)

  def draw_self_orbit(self, ax, E_data):
    R1 = E_data[2] + self.height
    T1 = self.time * self.velocity
    N = 1000.0

    t = [T1 * i / N for i in np.arange(0, N, 1)]

    X = np.array([R1 * np.cos(2 * np.pi * w / T1) for w in t])
    Y = np.array([R1 * np.sin(2 * np.pi * w / T1) for w in t])
    Z = np.array([w * 0 for w in t])

    ax.plot(X, Y, Z, label='Satellite-orbit', color='r')
    ax.legend()
    
  
  def return_data(self):
    return [self.mass, self.height, self.velocity, self.time, self.init_pos]




fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#Satellite_orbit = 0

E = Earth(ax)
S = Satellite(ax, E.return_data())

#print(S.return_data())
#print(E.return_data())

plt.show()