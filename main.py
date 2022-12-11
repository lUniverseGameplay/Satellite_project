import PIL
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from scipy.integrate import odeint
import math

def to_coord(c1, c2, rad):
    f_c1 = c1 * round(math.cos(math.radians(rad)), 5) - c2 * math.sin(math.radians(rad))
    f_c2 = c1 * math.sin(math.radians(rad)) + c2 * round(math.cos(math.radians(rad)), 5)
    return [f_c1, f_c2]

def odefun(lst: np.ndarray, t: float) -> np.ndarray:
    r = np.array([-1 * lst[0], -1 * lst[1], -1 * lst[2]])
    norm_r = np.linalg.norm(r) ** 3
    ax = G * M * -1 * lst[0] / norm_r + S.mass * -1 * lst[0] / norm_r
    ay = G * M * -1 * lst[1] / norm_r + S.mass * -1 * lst[1] / norm_r
    az = G * M * -1 * lst[2] / norm_r + S.mass * -1 * lst[2] / norm_r
    return np.array([lst[3], lst[4], lst[5], ax, ay, az])


class Earth():
  def __init__(self):
    global G, M, E_radius
    G = 6.673 * (10 ** (-11))
    M = 5.972 * (10 ** 24)
    E_radius = 6371000

    self.stratosphere_radius = E_radius + 50000
  
  def draw_me(self, ax):
    # Open Image in PIL
    bm = PIL.Image.open('data/earthicefreesm.jpg') 
    bm = np.array(bm.resize([d * 2 // 1 for d in bm.size]))/256

    # Create Earth with texture3
    lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180 
    lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180

    x = E_radius * np.outer(np.cos(lons), np.cos(lats)).T
    y = E_radius * np.outer(np.sin(lons), np.cos(lats)).T
    z = E_radius * np.outer(np.ones(np.size(lons)), np.sin(lats)).T

    # Plot the Earth
    ax.plot_surface(x, y, z, rstride=4, cstride=4, alpha=0.4, facecolors=bm)
    #ax.plot_surface(x, y, z, linewidth=1, alpha=0.3, color='g')
  

  def draw_stratosphere(self, ax):
    # Create Erath's stratosphere
    u_stratosphere = np.linspace(0, 2 * np.pi, 100)
    v_stratosphere = np.linspace(0, np.pi, 100)

    x_stratosphere = self.stratosphere_radius * np.outer(np.cos(u_stratosphere), np.sin(v_stratosphere))
    y_stratosphere = self.stratosphere_radius * np.outer(np.sin(u_stratosphere), np.sin(v_stratosphere))
    z_stratosphere = self.stratosphere_radius * np.outer(np.ones(np.size(u_stratosphere)), np.cos(v_stratosphere))

    # Plot the stratosphere
    ax.plot_wireframe(x_stratosphere, y_stratosphere, z_stratosphere, linewidth=1, alpha=0.3)
  
  def check_point_in_atmosphere(self, x, y, z):
    return S.height <= 50000

  def return_data(self):
    return [G, M, E_radius, self.stratosphere_radius, self.Oz]


class Satellite():
  def __init__(self):
    self.mass = 420000
    self.height = 437 * (10 ** 3)
    self.velocity = 7654
    #self.velocity = 7654
    self.time = 90 * 60

  def draw_self_orbit(self, ax, trajectory_corrected):
    X = np.array([i[0] for i in trajectory_corrected])
    Y = np.array([i[1] for i in trajectory_corrected])
    Z = np.array([i[2] for i in trajectory_corrected])

    ax.plot(X, Y, Z, label=str("The satellite is located in the Earth's atmosphere: "+str(E.check_point_in_atmosphere(x, y, z))), color='r', linewidth=1.5)
    ax.legend()
    
  
  def return_data(self):
    return [self.mass, self.height, self.velocity, self.time]




fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
#ax = fig.add_subplot(2, 2, 1, projection='3d')

#Satellite_orbit = 0

E = Earth()
S = Satellite()


#North = -43.07
#East = -61.5

#North = float(input())
North = -42.52

#East = float(input())
East = -10

x, z = to_coord(E_radius + S.height, 0, East)
x, y = to_coord(x, 0, North)

ax.scatter(int(x), int(y), int(z), marker='o', color='k')

coordinats = [int(x), int(y), int(z)]

tspan = np.linspace(0, 3 * S.time, 10 ** 5)

vel_list = [0, 0, S.velocity]

c_and_v = np.array(coordinats + vel_list)

xd = odeint(odefun, c_and_v, tspan)

trajectory = [i[0:3] for i in xd]
velocity = [i[3:6] for i in xd]

trajectory_corrected = np.zeros(np.shape(trajectory))

kinetic_enegry = np.zeros(np.shape(trajectory)[0])
potential_enegry = np.zeros(np.shape(trajectory)[0])

for i in range(len(tspan)):
  current_time = tspan[i]

  angle_Erth_rotation = -2 * np.pi * current_time / (24 * 60 * 60)

  current_point = np.array(trajectory[i][:])

  trajectory_corrected[i] = current_point

  kinetic_enegry[i] = 0.5 * S.mass * np.dot(velocity[i], velocity[0])
  potential_enegry[i] = -1 * G * M * S.mass / np.linalg.norm(current_point)

total_energy = potential_enegry + kinetic_enegry

E.draw_me(ax)
E.draw_stratosphere(ax)

S.draw_self_orbit(ax, trajectory_corrected)

#ax.scatter(0, 0, 0, marker='o', color='g')

ax.scatter(10 ** 7, 10 ** 7, 10 ** 7, marker='o', color='k', alpha=0)
ax.scatter(-10 ** 7, 10 ** 7, 10 ** 7, marker='o', color='k', alpha=0)
ax.scatter(10 ** 7, -10 ** 7, 10 ** 7, marker='o', color='k', alpha=0)
ax.scatter(-10 ** 7, -10 ** 7, 10 ** 7, marker='o', color='k', alpha=0)

#ax2 = fig.add_subplot(2, 2, 2)

#ax2.plot(tspan, kinetic_enegry, 'r', label="kinetic")
#ax2.plot(tspan, potential_enegry, 'b', label="potential")
#ax2.plot(tspan, total_energy, 'k', label="total")
#ax2.legend()
#ax2.set_title(['change in total energy: ', total_energy[-1] - total_energy[0]])

plt.show()