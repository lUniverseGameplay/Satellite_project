import PIL
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


class Earth():
  def __init__(self):
    self.G = 6.673 * (10 ** (-11))
    self.M = 5.972 * (10 ** 24)
    self.radius = 6371000
    self.Oz = [0, 0, 1]

    self.stratosphere_radius = self.radius + 50000
    
    self.draw_me()
    self.draw_stratosphere()
  
  def draw_me(self):
    # Open Image in PIL
    bm = PIL.Image.open('data/earthicefreesm.jpg') 
    bm = np.array(bm.resize([d * 2 // 1 for d in bm.size]))/256

    # Create Earth with texture
    lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180 
    lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180 


    fig = plt.figure()
    self.ax = fig.add_subplot(111, projection='3d')

    x = self.radius * np.outer(np.cos(lons), np.cos(lats)).T
    y = self.radius * np.outer(np.sin(lons), np.cos(lats)).T
    z = self.radius * np.outer(np.ones(np.size(lons)), np.sin(lats)).T

    # Plot the Earth
    self.ax.plot_surface(x, y, z, rstride=4, cstride=4, facecolors = bm)
  

  def draw_stratosphere(self):
    # Create Erath's stratosphere
    u_stratosphere = np.linspace(0, 2 * np.pi, 100)
    v_stratosphere = np.linspace(0, np.pi, 100)

    x_stratosphere = self.stratosphere_radius * np.outer(np.cos(u_stratosphere), np.sin(v_stratosphere))
    y_stratosphere = self.stratosphere_radius * np.outer(np.sin(u_stratosphere), np.sin(v_stratosphere))
    z_stratosphere = self.stratosphere_radius * np.outer(np.ones(np.size(u_stratosphere)), np.cos(v_stratosphere))

    # Plot the stratosphere
    self.ax.plot_wireframe(x_stratosphere, y_stratosphere, z_stratosphere, linewidth=1, alpha=0.3)

  def return_data(self):
    return self.G, self.M, self.radius, self.stratosphere_radius, self.Oz


class Satellite():
  def __init__(self):
    self.mass = 457 * (10 ** 3)
    self.height = 0
    self.velocity = 7654
    self.time = 90 * 60
    self.North = float(input()) * np.pi / 180
    self.East = float(input())  * np.pi / 180

    self.init_pos = [np.cos(self.North) * np.cos(self.East),
    np.cos(self.North) * np.sin(self.East),
    np.sin(self.North)]

  def return_data(self):
    return self.mass, self.height, self.velocity, self.time, self.init_pos




#Satellite_orbit = 0

S = Satellite()
E = Earth()

print(S.return_data())
print(E.return_data())

plt.show()