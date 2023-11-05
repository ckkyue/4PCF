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

def A(l, m, primary_vert, secondary_vert, bin, weight_lists, avg_cal=True, config_space_avg=False):
    # Get the number of vertices
    num_vert = secondary_vert.shape[0]
    # Initialize sum variable
    sum = 0
    # Calculate relative positions in spherical coordinates
    rel_pos_sphe = cart_to_sphe(secondary_vert - primary_vert)
    # Check conditions for calculation
    if avg_cal and not config_space_avg:
        # Iterate over vertices
        for i in range(num_vert):
            # Check if the radius is within the specified bin range
            if bin[0] < rel_pos_sphe[i][0] < bin[1]:
                # Accumulate the weighted sum using spherical harmonic function
                sum += weight_lists[i] * spe.sph_harm(m, l, rel_pos_sphe[i][2], rel_pos_sphe[i][1])
        # Calculate and return the average
        return sum / shell_vol(bin[0], bin[1])

def estimator(l1, l2, l3, Data_catalog, Data_catalog_weight_lists, bin_lists, Random_catalog, Random_catalog_weight_lists, DmR_status, space_volume, cache_flag=0, avg_cal=True, config_space_avg=False):
    sum = 0
    # DDDD
    if DmR_status == 0:
        num_data = Data_catalog.shape[0]
        num_weights = Data_catalog_weight_lists.shape[0]
        assert num_data == num_weights

    for i in range(num_data):
        primary = Data_catalog[i]
        secondary = np.delete(Data_catalog, i, axis=0)
        weight_lists_primary = Data_catalog_weight_lists[i]
        weight_lists_secondary = np.delete(Data_catalog_weight_lists, i, axis=0)

        for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
            sum += weight_lists_primary * C(l1, l2, l3, m1, m2, m3) * A(l1, m1, primary, secondary, bin_lists[0], weight_lists_secondary, avg_cal, config_space_avg) \
                   * A(l2, m2, primary, secondary, bin_lists[0], weight_lists_secondary, avg_cal, config_space_avg) \
                   * A(l3, m3, primary, secondary, bin_lists[0], weight_lists_secondary, avg_cal, config_space_avg)

# End of the program
end = datetime.now()
print(end-start)