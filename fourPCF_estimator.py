import numpy as np
from sympy.physics.wigner import wigner_3j
import scipy.special as spe
import math
from datetime import datetime

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

def C(l1, l2, l3, m1, m2, m3):
    return (-1)**(l1+l2+l3)*wigner_3j(l1, l2, l3, m1, m2, m3)

def a(l, m, primary_vert, secondary_vert, bins, weights, avg_cal=True, config_space_avg=False):
    # Get the number of vertices
    num_vert = secondary_vert.shape[0]
    # Initialize sum variable
    sum = 0
    # Calculate relative positions in spherical coordinates
    rel_pos_sphe = cart_to_sphe(secondary_vert - primary_vert)

    if (avg_cal == True) and (config_space_avg == True):
        for i in range(num_vert):
            if bins[0] < rel_pos_sphe[i][0] < bins[1]:
                sum += weights[i] * spe.sph_harm(m, l, rel_pos_sphe[i][2], rel_pos_sphe[i][1])
        return sum / shell_vol(bins[0], bins[1])
            
    elif (avg_cal == True) and (config_space_avg == False):
        for i in range(num_vert):
            if bins[0] < rel_pos_sphe[i][0] < bins[1]:
                sum += rel_pos_sphe[i][0] * weights[i] * spe.sph_harm(m, l, rel_pos_sphe[i][2], rel_pos_sphe[i][1])
        return sum / shell_vol(bins[0], bins[1])
            
    elif (avg_cal == False) and (config_space_avg == False):
        for i in range(num_vert):
            if bins[0] < rel_pos_sphe[i][0] < bins[1]:
                sum += rel_pos_sphe[i][0] * weights[i] * spe.sph_harm(m, l, rel_pos_sphe[i][2], rel_pos_sphe[i][1])
        return sum
            
    else:        
        for i in range(num_vert):
            if bins[0] < rel_pos_sphe[i][0] < bins[1]:
                sum += weights[i] * spe.sph_harm(m, l, rel_pos_sphe[i][2], rel_pos_sphe[i][1])
        return sum
    
def estimator(l1, l2, l3, Data_catalog, Data_catalog_weights, bins, Random_catalog, Random_catalog_weights, DmR_status, space_vol, cache_flag=0, avg_cal=True, config_space_avg=False):
    sum = 0

    # DDDD
    if DmR_status == 0:
        num_data = Data_catalog.shape[0]
        num_weights = Data_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Data_catalog[i]
            secondary = np.delete(Data_catalog, i, axis=0)
            weights_primary = Data_catalog_weights[i]
            weights_secondary = np.delete(Data_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, secondary, bins[0], weights_secondary, avg_cal, config_space_avg) \
                * a(l2, m2, primary, secondary, bins[1], weights_secondary, avg_cal, config_space_avg) \
                * a(l3, m3, primary, secondary, bins[2], weights_secondary, avg_cal, config_space_avg)
    
    # DDDR
    elif DmR_status == 1:
        num_data = Data_catalog.shape[0]
        num_weights = Data_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Data_catalog[i]
            secondary = np.delete(Data_catalog, i, axis=0)
            weights_primary = Data_catalog_weights[i]
            weights_secondary = np.delete(Data_catalog_weights, i, axis=0)
            
            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, secondary, bins[0], weights_secondary) \
                * a(l2, m2, primary, secondary, bins[1], weights_secondary) \
                * a(l3, m3, primary, Random_catalog, bins[2], Random_catalog_weights)

    # DDRD
    elif DmR_status == 2:
        num_data = Data_catalog.shape[0]
        num_weights = Data_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Data_catalog[i]
            secondary = np.delete(Data_catalog, i, axis=0)
            weights_primary = Data_catalog_weights[i]
            weights_secondary = np.delete(Data_catalog_weights, i, axis=0)
            
            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, secondary, bins[0], weights_secondary) \
                * a(l2, m2, primary, Random_catalog, bins[1], Random_catalog_weights) \
                * a(l3, m3, primary, secondary, bins[2], weights_secondary)
                
    # DRDD
    elif DmR_status == 3:
        num_data = Data_catalog.shape[0]
        num_weights = Data_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Data_catalog[i]
            secondary = np.delete(Data_catalog, i, axis=0)
            weights_primary = Data_catalog_weights[i]
            weights_secondary = np.delete(Data_catalog_weights, i, axis=0)
            
            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, Random_catalog, bins[0], Random_catalog_weights) \
                * a(l2, m2, primary, secondary, bins[1], weights_secondary) \
                * a(l3, m3, primary, secondary, bins[2], weights_secondary)

    # RDDD
    elif DmR_status == 4:
        num_data = Random_catalog.shape[0]
        num_weights = Random_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Random_catalog[i]
            secondary = np.delete(Random_catalog, i, axis=0)
            weights_primary = Random_catalog_weights[i]
            weights_secondary = np.delete(Random_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, secondary, bins[0], Data_catalog_weights) \
                * a(l2, m2, primary, secondary, bins[1], Data_catalog_weights) \
                * a(l3, m3, primary, secondary, bins[2], Data_catalog_weights)

    # DDRR
    elif DmR_status == 5:
        num_data = Data_catalog.shape[0]
        num_weights = Data_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Data_catalog[i]
            secondary = np.delete(Data_catalog, i, axis=0)
            weights_primary = Data_catalog_weights[i]
            weights_secondary = np.delete(Data_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, secondary, bins[0], weights_secondary) \
                * a(l2, m2, primary, Random_catalog, bins[1], Random_catalog_weights) \
                * a(l3, m3, primary, Random_catalog, bins[2], Random_catalog_weights)

    # DRRD
    elif DmR_status == 6:
        num_data = Data_catalog.shape[0]
        num_weights = Data_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Data_catalog[i]
            secondary = np.delete(Data_catalog, i, axis=0)
            weights_primary = Data_catalog_weights[i]
            weights_secondary = np.delete(Data_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, Random_catalog, bins[0], Random_catalog_weights) \
                * a(l2, m2, primary, Random_catalog, bins[1], Random_catalog_weights) \
                * a(l3, m3, primary, secondary, bins[2], weights_secondary)

    # DRDR
    elif DmR_status == 7:
        num_data = Data_catalog.shape[0]
        num_weights = Data_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Data_catalog[i]
            secondary = np.delete(Data_catalog, i, axis=0)
            weights_primary = Data_catalog_weights[i]
            weights_secondary = np.delete(Data_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, Random_catalog, bins[0], Random_catalog_weights) \
                * a(l2, m2, primary, secondary, bins[1], weights_secondary) \
                * a(l3, m3, primary, Random_catalog, bins[2], Random_catalog_weights)

    # RRDD
    elif DmR_status == 8:
        num_data = Random_catalog.shape[0]
        num_weights = Random_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Random_catalog[i]
            secondary = np.delete(Random_catalog, i, axis=0)
            weights_primary = Random_catalog_weights[i]
            weights_secondary = np.delete(Random_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, Random_catalog, bins[0], Random_catalog_weights) \
                * a(l2, m2, primary, Data_catalog, bins[1], Data_catalog_weights) \
                * a(l3, m3, primary, Data_catalog, bins[2], Data_catalog_weights)

    # RDRD
    elif DmR_status == 9:
        num_data = Random_catalog.shape[0]
        num_weights = Random_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Random_catalog[i]
            secondary = np.delete(Random_catalog, i, axis=0)
            weights_primary = Random_catalog_weights[i]
            weights_secondary = np.delete(Random_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, Data_catalog, bins[0], Data_catalog_weights) \
                * a(l2, m2, primary, secondary, bins[1], weights_secondary) \
                * a(l3, m3, primary, Data_catalog, bins[2], Data_catalog_weights)

    # RDDR
    elif DmR_status == 10:
        num_data = Random_catalog.shape[0]
        num_weights = Random_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Random_catalog[i]
            secondary = np.delete(Random_catalog, i, axis=0)
            weights_primary = Random_catalog_weights[i]
            weights_secondary = np.delete(Random_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, Data_catalog, bins[0], Data_catalog_weights) \
                * a(l2, m2, primary, Data_catalog, bins[1], Data_catalog_weights) \
                * a(l3, m3, primary, secondary, bins[2], weights_secondary)

    # DRRR
    elif DmR_status == 11:
        num_data = Data_catalog.shape[0]
        num_weights = Data_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Data_catalog[i]
            secondary = np.delete(Data_catalog, i, axis=0)
            weights_primary = Data_catalog_weights[i]
            weights_secondary = np.delete(Data_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, Random_catalog, bins[0], Random_catalog_weights) \
                * a(l2, m2, primary, Random_catalog, bins[1], Random_catalog_weights) \
                * a(l3, m3, primary, Random_catalog, bins[2], Random_catalog_weights)

    # RDRR
    elif DmR_status == 12:
        num_data = Random_catalog.shape[0]
        num_weights = Random_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Random_catalog[i]
            secondary = np.delete(Random_catalog, i, axis=0)
            weights_primary = Random_catalog_weights[i]
            weights_secondary = np.delete(Random_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, Data_catalog, bins[0], Data_catalog_weights) \
                * a(l2, m2, primary, secondary, bins[1], weights_secondary) \
                * a(l3, m3, primary, secondary, bins[2], weights_secondary)

    # RRDR
    elif DmR_status == 13:
        num_data = Random_catalog.shape[0]
        num_weights = Random_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Random_catalog[i]
            secondary = np.delete(Random_catalog, i, axis=0)
            weights_primary = Random_catalog_weights[i]
            weights_secondary = np.delete(Random_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, secondary, bins[0], weights_secondary) \
                * a(l2, m2, primary, Data_catalog, bins[1], Data_catalog_weights) \
                * a(l3, m3, primary, secondary, bins[2], weights_secondary)

    # RRRD
    elif DmR_status == 14:
        num_data = Random_catalog.shape[0]
        num_weights = Random_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Random_catalog[i]
            secondary = np.delete(Random_catalog, i, axis=0)
            weights_primary = Random_catalog_weights[i]
            weights_secondary = np.delete(Random_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, secondary, bins[0], weights_secondary) \
                * a(l2, m2, primary, secondary, bins[1], weights_secondary) \
                * a(l3, m3, primary, Data_catalog, bins[2], Data_catalog_weights)

    # RRRR
    elif DmR_status == 15:
        num_data = Random_catalog.shape[0]
        num_weights = Random_catalog_weights.shape[0]
        assert num_data == num_weights

        for i in range(num_data):
            primary = Random_catalog[i]
            secondary = np.delete(Random_catalog, i, axis=0)
            weights_primary = Random_catalog_weights[i]
            weights_secondary = np.delete(Random_catalog_weights, i, axis=0)

            for m1, m2, m3 in zip(range(-l1, l1+1), range(-l2, l2+1), range(-l3, l3+1)):
                sum += weights_primary * C(l1, l2, l3, m1, m2, m3) \
                * a(l1, m1, primary, secondary, bins[0], weights_secondary) \
                * a(l2, m2, primary, secondary, bins[1], weights_secondary) \
                * a(l3, m3, primary, secondary, bins[2], weights_secondary)
    return sum/space_vol