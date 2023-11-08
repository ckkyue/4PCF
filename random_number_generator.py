import numpy as np
import random

def min_sep_test(point, points, min_sep):
    result = True
    for i in range(len(points)):
        if np.linalg.norm(point - points[i]) < min_sep:
            result = False
            break
    return result

def generate_1d_random(num, space, min_sep):
    points = np.empty((0, 1))
    x_min = space[0][0]
    x_max = space[0][1]

    while num > 0:
        x = random.uniform(x_min, x_max)

        point = np.array([x])
        if min_sep_test(point, points, min_sep):
            points = np.append(points, point.reshape((1, 1)), axis=0)
            num -= 1
    
    return points

def generate_3d_random(num, space, min_sep):
    points = np.empty((0, 3))
    x_min = space[0][0]
    x_max = space[0][1]
    y_min = space[1][0]
    y_max = space[1][1]
    z_min = space[2][0]
    z_max = space[2][1]

    while num > 0:
        x = random.uniform(x_min, x_max)
        y = random.uniform(y_min, y_max)
        z = random.uniform(z_min, z_max)

        point = np.array([x, y, z])
        if min_sep_test(point, points, min_sep):
            points = np.append(points, point.reshape((1, 3)), axis=0)
            num -= 1

    return points