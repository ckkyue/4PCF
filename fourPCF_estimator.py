import numpy as np
import math
from read_data import read, get_weights, get_carts
from sympy.physics.wigner import wigner_3j
import scipy.special as spe
from tqdm import tqdm

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

# Volume of spherical shell
def shell_vol(bin_min, bin_max):
    return 4/3*np.pi*(bin_max**3 - bin_min**3)

def a(l, m, primary_vert, secondary_vert, bin_min, bin_max, weight):
    # Get the number of vertices
    num_vert = secondary_vert.shape[0]
    # Initialize sum variable
    sum = 0
    v_b = 1
    # Calculate relative positions in spherical coordinates
    rel_pos_sphe = cart_to_sphe(secondary_vert - primary_vert)
    for i in range(num_vert):
        if bin_min < rel_pos_sphe[i][0] < bin_max:
            sum += weight * spe.sph_harm(m, l, rel_pos_sphe[i][2], rel_pos_sphe[i][1])
    v_b = v_b * shell_vol(bin_min, bin_max)
    return sum / v_b

def estimator(l1, l2, l3, vertices, bins_min, bins_max, weights):
    num_vert = vertices.shape[0]
    sum = 0
    for i in range(num_vert):
        weight = weights[i]
        primary = vertices[i]
        secondary = np.delete(vertices, i, axis=0)
        for m1 in range(-l1, l1+1):
            for m2 in range(-l2, l2+1):
                for m3 in range(-l3, l3+1):
                    sum += np.float64(wigner_3j(l1, l2, l3, m1, m2, m3)) \
                        * a(l1, m1, primary, secondary, bins_min[0], bins_max[0], weight) \
                        * a(l2, m2, primary, secondary, bins_min[1], bins_max[1], weight) \
                        * a(l3, m3, primary, secondary, bins_min[2], bins_max[2], weight)
    return sum

# Extract data
# data_path = "test_sample.txt"
# data = read(data_path)
# vertices = get_carts(data)
# np.save("vertices_sample.npy", vertices)
# weights = get_weights(data)
# np.save("weights_sample.npy", weights)

vertices = np.load("vertices_sample.npy")
weights = np.load("weights_sample.npy")

radial_bins = []
for i in range(10):
    r1 = np.linspace(20, 160, 10 + 1)[i]
    for j in range(i + 1, 10):
        r2 = np.linspace(20, 160, 10 + 1)[j]
        if r2 - r1 > 14:
            for k in range(j + 1, 10):
                r3 = np.linspace(20, 160, 10 + 1)[k]
                if r3 - r2 > 14:
                    radial_bins.append([r1, r2, r3])

# Estimate 4PCF of test sample
zeta = []
for bins in tqdm(radial_bins, desc="Processing bins"):
    bins_min = np.array(bins)
    bins_max = bins_min + 14
    zeta_bin = np.imag(estimator(1, 1, 1, vertices, bins_min, bins_max, weights))
    zeta.append(zeta_bin)
np.save("zeta_test_sample.npy", zeta)
print(zeta)