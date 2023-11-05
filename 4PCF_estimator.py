import numpy as np
import wigner as wg
import scipy.special as spe
import math
from datetime import datetime

# Start of the program
start = datetime.now()

# Longitude of a point within [0, 2*pi) given two 1D arrays
def longitude(x, y):
    output = np.empty(len(x))
    for i in range(len(x)):
        value = math.atan2(y[i], x[i])
        if value < 0:
            value += 2 * np.pi
        output[i] = value  
    return output

# Convert cartesian to spherical coordinates
def cart_to_sphe(cart):
    x, y, z = cart[:, 0], cart[:, 1], cart[:, 2]
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    spherical_coords = np.array([r, theta, phi]).T
    spherical_coords[np.isnan(spherical_coords)] = 0
    return spherical_coords

# Bin function
def bin_func(x, bin_min, bin_max):
    return int(bin_min < x < bin_max)

# Volume of sphere
def sphe_vol(r):
    return 4/3*np.pi*r**3

# Volume of spherical shell
def shell_vol(bin_min, bin_max):
    return sphe_vol(bin_min) - sphe_vol(bin_max)

# wigner_3j symbol
def wigner_3j(l1, l2, l3, m1, m2, m3):
    if m1 + m2 + m3 != 0:
        return 0
    intermedia = wg.wigner_3jm(l1, l2, l3, m1)
    index = int(m2 - intermedia[0])
    return intermedia[2][index]

def C(l1, l2, l3, m1, m2, m3):
    return (-1)*(l1+l2+l3)*wigner_3j(l1, l2, l3, m1, m2, m3)

def A_func(l, m, primary_vert, secondary_vert, bin, weight_lists, avg_cal = True, config_space_avg = False):
    num_vert = secondary_vert.shape[0]
    sum = 0
    b_min = bin[0]
    b_max = bin[1]
    rel_pos = secondary_vert - primary_vert
    rel_pos_sphe = cart_to_sphe(rel_pos)

    if (avg_cal == True) and (config_space_avg == False):
        for i in range(num_vert):
            if rel_pos_sphe[i][0] < b_max and rel_pos_sphe[i][0] > b_min:
                sum += rel_pos_sphe[i][0] * spe.sph_harm(m, l, rel_pos_sphe[i][2], rel_pos_sphe[i][1])
        return sum / shell_vol(b_min, b_max)

# End of the program
end = datetime.now()
print(end-start)