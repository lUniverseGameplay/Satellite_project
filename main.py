import PIL
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import math


def odefun(x: np.ndarray, t: float):
  p = np.array([np.array(x[3:6]), np.array([-1 * 6.673 * (10 ** (-11)) * 5.972 * (10 ** 24) * i / (round((np.linalg.norm(x[:3])), 4) ** 3) for i in x[0:3]])])
  return [p[0][0], p[0][1], p[0][2], p[1][0], p[1][1], p[1][2]]
  
def rotz(gamma):
  return [[np.cos(gamma), np.sin(gamma) * -1, 0], [np.sin(gamma), np.cos(gamma), 0], [0, 0, 1]]

def get_orbit(r, sol):
    phi = 51.6 * np.pi / 180
    p1 = -r[1] / r[0]
    p2 = -np.cos(phi) * r[2] / r[0]
    a = p1 ** 2 + 1
    b = 2 * p1 * p2
    c = p2 ** 2 - np.sin(phi) * np.sin(phi)
    y1 = (-b + math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    y2 = (-b - math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    x1 = p1 * y1 + p2
    x2 = p1 * y2 + p2
    z = np.cos(phi)
    n1 = [x1, y1, z]
    n2 = [x2, y2, z]
    if sol:
        return n1
    else:
        return n2


class Earth():
  def __init__(self, ax):
    global G, M, radius
    G = 6.673 * (10 ** (-11))
    M = 5.972 * (10 ** 24)
    radius = 6371000
    self.Oz = [0, 0, 1]

    self.stratosphere_radius = radius + 50000
    
    #self.draw_me(ax)
    #self.draw_stratosphere(ax)
  
  def draw_me(self, ax):
    # Open Image in PIL
    bm = PIL.Image.open('data/earthicefreesm.jpg') 
    bm = np.array(bm.resize([d * 2 // 1 for d in bm.size]))/256

    # Create Earth with texture3
    lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180 
    lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180

    x = radius * np.outer(np.cos(lons), np.cos(lats)).T
    y = radius * np.outer(np.sin(lons), np.cos(lats)).T
    z = radius * np.outer(np.ones(np.size(lons)), np.sin(lats)).T

    # Plot the Earth
    ax.plot_surface(x, y, z, rstride=4, cstride=4, facecolors=bm)
  

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
    return [G, M, radius, self.stratosphere_radius, self.Oz]


class Satellite():
  def __init__(self, ax):
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

    r0 = np.asarray([i * (radius + self.height) for i in self.init_pos])
    v0 = self.tau * self.velocity

    tspan = np.linspace(0, 3 * self.time, 10 ** 5)
    x0 = [int(r0[0]), int(r0[1]), int(r0[2]), int(v0[0]), int(v0[1]), int(v0[2])]

    #odefun = lambda t, x: [np.asarray(x[1]), np.asarray([-1 * E_data[0] * E_data[1] * i / ((np.linalg.norm(i)) ** 3) for i in x[0]])]

    #xd = odeint(odefun, x0, tspan)

    self.x = odeint(odefun, x0, tspan)

    trajectory = [i[0:3] for i in self.x]
    velocity = [i[3:6] for i in self.x]

    trajectory_corrected = np.zeros(np.shape(trajectory))

    kinetic_enegry = np.zeros(np.shape(trajectory)[0])
    potential_enegry = np.zeros(np.shape(trajectory)[0])

    file_log = open('data/trajectory_coorected.txt', 'w')

    for i in range(len(tspan)):
      current_time = tspan[i]
      #print(current_time)
      angle_Erth_rotation = -2 * np.pi * current_time / (24 * 3600)
      #print(angle_Erth_rotation)

      current_point = np.array(trajectory[i][:]).transpose()
      current_point_corrected = rotz(angle_Erth_rotation) * current_point

      trajectory_corrected[i] = [current_point_corrected[i][i] for i in range(len(current_point_corrected))]

      file_log.write(str([int(self.x[i][0]), int(self.x[i][1]), int(self.x[i][2]), int(self.x[i][3]), int(self.x[i][4]), int(self.x[i][5])]))
      file_log.write('\n')

      kinetic_enegry[i] = 0.5 * self.mass * np.dot(velocity[i], velocity[0])
      potential_enegry[i] = -1 * G * M * self.mass / np.linalg.norm(current_point)
      #if i == 2493 or i == 2492 or i == 2494:
        #print(i, current_time, angle_Erth_rotation, current_point, current_point_corrected, trajectory_corrected[i], sep='\n')
        #print('---------------------------------------------------------------------')
    
    file_log.close()
    
    total_energy = potential_enegry + kinetic_enegry

    self.draw_self_orbit(ax, trajectory_corrected)

  def draw_self_orbit(self, ax, trajectory_corrected):
    X = np.array([i[0] for i in trajectory_corrected])
    Y = np.array([i[1] for i in trajectory_corrected])
    Z = np.array([i[2] for i in trajectory_corrected])

    ax.plot(X, Y, Z, label='Satellite-orbit', color='r', linewidth=1.5)
    ax.scatter(int(self.x[0][0]), int(self.x[0][1]), int(self.x[0][2]), marker='o', color='k')
    ax.legend()
    
  
  def return_data(self):
    return [self.mass, self.height, self.velocity, self.time, self.init_pos]




fig = plt.figure()
ax = fig.add_subplot(projection='3d')

#Satellite_orbit = 0

E = Earth(ax)
S = Satellite(ax)

#print(S.return_data())
#print(E.return_data())

plt.show()