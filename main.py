import PIL
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from scipy.integrate import odeint
import math


def to_coord(c1, c2, rad):
  f_c1 = c1 * math.cos(math.radians(rad)) - c2 * math.sin(math.radians(rad))
  f_c2 = c1 * math.sin(math.radians(rad)) + c2 * math.cos(math.radians(rad))
  return [f_c1, f_c2]

def odefun(lst: np.ndarray, t: float) -> np.ndarray:
  r = np.array([-1 * lst[0], -1 * lst[1], 1 * lst[2]])
  norm_r = np.linalg.norm(r) ** 3
  ax = G * (M * -1 * lst[0] / norm_r + S.mass * -1 * lst[0] / norm_r)
  ay = G * (M * -1 * lst[1] / norm_r + S.mass * -1 * lst[1] / norm_r)
  az = G * (M * -1 * lst[2] / norm_r + S.mass * 1 * lst[2] / norm_r)
  return np.array([lst[3], lst[4], lst[5], ax, ay, az])


class Earth():
  def __init__(self):
    global G, M, E_radius
    G = 6.673 * (10 ** (-11))
    M = 5.972 * (10 ** 24)
    E_radius = 6371000

    self.stratosphere_high = 50000
    self.orbit_min = 160000
    self.orbit_max = 2 * 10 ** 6
  
  def draw_me(self, ax, East):
    # Open Image in PIL
    bm = PIL.Image.open('data/earthicefreesm_true_size.jpg') 
    bm = np.array(bm.resize([d * 2 // 1 for d in bm.size]))/256

    # Create Earth with texture3
    lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180 
    lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180

    x = E_radius * np.outer(np.cos(lons), np.cos(lats)).T
    y = E_radius * np.outer(np.sin(lons), np.cos(lats)).T
    z = E_radius * np.outer(np.ones(np.size(lons)), np.sin(lats)).T

    a = np.pi / (180 / -East)
    t = np.transpose(np.array([x,y,z]), (1,2,0))
    m = [[np.cos(a), 0, np.sin(a)],[0,1,0],[-np.sin(a), 0, np.cos(a)]]
    x,y,z = np.transpose(t @ m, (2,0,1))

    ax.plot_surface(x, y, z, rstride=4, cstride=4, alpha=0.4, facecolors=bm)
  

  def draw_stratosphere(self, ax):
    # Create Erath's stratosphere
    u_stratosphere = np.linspace(0, 2 * np.pi, 100)
    v_stratosphere = np.linspace(0, np.pi, 100)

    x_stratosphere = (self.stratosphere_high + E_radius) * np.outer(np.cos(u_stratosphere), np.sin(v_stratosphere))
    y_stratosphere = (self.stratosphere_high + E_radius) * np.outer(np.sin(u_stratosphere), np.sin(v_stratosphere))
    z_stratosphere = (self.stratosphere_high + E_radius) * np.outer(np.ones(np.size(u_stratosphere)), np.cos(v_stratosphere))

    # Plot the stratosphere
    ax.plot_wireframe(x_stratosphere, y_stratosphere, z_stratosphere, linewidth=1, alpha=0.3)
  
  def check_point_in_atmosphere(self, x, y, z):
    f1 = (x ** 2 + y ** 2 + 10 ** 13) <= (self.stratosphere_high + E_radius) ** 2 # ???????????????? ???? ?????????????????????? ?? ????????????????????????
    f2 = (x ** 2 + z ** 2 + 10 ** 13) <= (self.stratosphere_high + E_radius) ** 2 # ???????????????? ???? ?????????????????????? ?? ????????????????????????
    f3 = (y ** 2 + z ** 2 + 10 ** 13) <= (self.stratosphere_high + E_radius) ** 2 # ???????????????? ???? ?????????????????????? ?? ????????????????????????
    return  f1 and f2 and f3
  
  def check_point_in_Earth(self, x, y, z):
    f1 = (x ** 2 + y ** 2 + 10 ** 13) <= (E_radius) ** 2 # ???????????????? ???? ?????????????????????? ?? ????????????????????????
    f2 = (x ** 2 + z ** 2 + 10 ** 13) <= (E_radius) ** 2 # ???????????????? ???? ?????????????????????? ?? ????????????????????????
    f3 = (y ** 2 + z ** 2 + 10 ** 13) <= (E_radius) ** 2 # ???????????????? ???? ?????????????????????? ?? ????????????????????????
    return  f1 and f2 and f3

  def check_point_out_orbit(self, x, y, z):
    f1 = (x ** 2 + y ** 2 + 10 ** 13) >= ((self.orbit_min + E_radius) ** 2) and (x ** 2 + y ** 2 + 10 ** 13) <= ((self.orbit_max + E_radius) ** 2) # ???????????????? ???? ?????????????????????? ?? ????????????????????????
    f2 = (x ** 2 + z ** 2 + 10 ** 13) >= ((self.orbit_min + E_radius) ** 2) and (x ** 2 + z ** 2 + 10 ** 13) <= ((self.orbit_max + E_radius) ** 2) # ???????????????? ???? ?????????????????????? ?? ????????????????????????
    f3 = (y ** 2 + z ** 2 + 10 ** 13) >= ((self.orbit_min + E_radius) ** 2) and (y ** 2 + z ** 2 + 10 ** 13) <= ((self.orbit_max + E_radius) ** 2) # ???????????????? ???? ?????????????????????? ?? ????????????????????????
    return  f1 or f2 or f3

  def return_data(self):
    return [G, M, E_radius, self.stratosphere_high, self.Oz]


class Satellite():
  def __init__(self):
    self.mass = 420000
    self.height = 437 * (10 ** 3)
    self.velocity = 8000
    self.velocity2 = 10000
    self.time = 90 * 60
    self.time2 = 90 * 60 * 3

    self.color_orbit = ['r', 'g']

    self.in_e_orbit = False
    self.in_e_stratosphere = False
    self.in_Earth = False

  def draw_self_orbit(self, ax, trajectory_corrected, ind):
    X = np.array([i[0] for i in trajectory_corrected])
    Y = np.array([i[1] for i in trajectory_corrected])
    Z = np.array([i[2] for i in trajectory_corrected])

    label_name = "The satellite is located in the Earth's atmosphere: " + str(self.in_e_stratosphere) + '\n'
    label_name += "The satellite is located in the Earth's orbit: " + str(self.in_e_orbit) + '\n'
    label_name += "The satellite has fallen: " + str(self.in_Earth)
    
    ax.plot(X, Y, Z, label=label_name, color=self.color_orbit[ind], linewidth=1.5)
    ax.legend()
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')
    
  
  def return_data(self):
    return [self.mass, self.height, self.velocity, self.time]




fig = plt.figure(figsize=plt.figaspect(0.5)*1.5)
ax = fig.add_subplot(2, 2, 1, projection='3d')

E = Earth()
S = Satellite()

North = float(input('North coords\n')) % 360

East = float(input('\nEast coords\n')) % 90

East_standart = 0

x_st, z_st = to_coord(E_radius + S.height, 0, East_standart)
x_st, y_st = to_coord(x_st, 0, North)

ax.scatter(int(x_st), int(y_st), int(z_st), marker='o', color='k')

E.draw_me(ax, East)
E.draw_stratosphere(ax)
  
coordinats = [int(x_st), int(y_st), int(z_st)]

for j in range(2):
  if j:
    tspan = np.linspace(0, 3 * S.time2, 10 ** 5)
    vel_list = [0, 0, S.velocity2]
  else:
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

    current_point = np.array(trajectory[i][:]).transpose()

    trajectory_corrected[i] = current_point
    if not S.in_e_stratosphere:
      if E.check_point_in_atmosphere(math.fabs(current_point[0]), math.fabs(current_point[1]), math.fabs(current_point[2])):
        S.in_e_stratosphere = True
    if not S.in_Earth:
      if E.check_point_in_Earth(math.fabs(current_point[0]), math.fabs(current_point[1]), math.fabs(current_point[2])):
        S.in_Earth = True
    if not S.in_e_orbit:
      if E.check_point_out_orbit(math.fabs(current_point[0]), math.fabs(current_point[1]), math.fabs(current_point[2])):
        S.in_e_orbit = True

    kinetic_enegry[i] = 0.5 * S.mass * np.dot(velocity[i][:], velocity[i][:])
    potential_enegry[i] = -1 * G * M * S.mass / np.linalg.norm(current_point)

  total_energy = potential_enegry + kinetic_enegry

  S.draw_self_orbit(ax, trajectory_corrected, j)

  ax2 = fig.add_subplot(2, 2, 2 + j)

  ax2.plot(tspan, kinetic_enegry, 'r', label="kinetic")
  ax2.plot(tspan, potential_enegry, 'b', label="potential")
  ax2.plot(tspan, total_energy, 'k', label="total")
  ax2.legend()
  ax2.set_title(['change in total energy: ' + S.color_orbit[j] + ' orbit', total_energy[-1] - total_energy[0]])

ax.axis('scaled')
plt.show()