import numpy as np
import math

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