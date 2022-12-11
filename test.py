from numpy import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def f(y, t):
         y1, y2, y3, y4, y5, y6 = y
         return [y2, -(4*pi*pi*y1)/(y1**2+y3**2 +y5**2)**(3/2),y4,-(4*pi*pi*y3)/(y1**2+y3**2 +y5**2)**(3/2),y6,-(4*pi*pi*y5)/(y1**2+y3**2 +y5**2)**(3/2)]


t = linspace(0,300,10001)
y0 = [0.325514,-9.096111, -0.459460,-6.916686,0.166229,-1.305721]
[y1,y2, y3, y4,y5,y6] = odeint(f, y0, t, full_output=False).T
fig, ax = plt.subplots()
plt.title("Орбита кометы Галлея(расстояние в а.е., время в годах) \n Солнце в центре координат")
plt.xlabel('x(t)')
plt.ylabel('y(t)')
fig.set_facecolor('white')
ax.plot(y1,y3,linewidth=1)
circle = Circle((0, 0), 0.2, facecolor='orange')   
ax.add_patch(circle)
plt.axis([1,-21,-1,29])
plt.grid(True)
fig, ax = plt.subplots()
plt.title("Орбита кометы Галлея \n Солнце в центре координат")
plt.xlabel('x(t)')
plt.ylabel('z(t)')
fig.set_facecolor('white')
ax.plot(y1,y5,linewidth=1)
circle = Circle((0, 0), 0.1, facecolor='orange')   
ax.add_patch(circle)
plt.axis([1,-21,1,-11])
plt.grid(True)
fig, ax = plt.subplots()
plt.title("Орбита кометы Галлея \n Солнце в центре координат")
plt.xlabel('y(t)')
plt.ylabel('z(t)')
fig.set_facecolor('white')
ax.plot(y3,y5,linewidth=1)
circle = Circle((0, 0), 0.2, facecolor='orange')   
ax.add_patch(circle)
plt.axis([-1,29,1,-11])
plt.grid(True)
fig, ax = plt.subplots()
plt.title("Проекция скорости движения  кометы Галлея \n на плоскости ZOX и ZOY ")
ax.plot(t,y1,linewidth=1)
ax.plot(t,y3,linewidth=1)
plt.show()