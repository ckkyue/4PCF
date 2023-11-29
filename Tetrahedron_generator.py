import numpy as np
import random
from random_number_generator import generate_1d_random, generate_3d_random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import fourPCF_estimator as fe
from tqdm import tqdm

folders = ["Figure"]
for folder in folders:
    if not os.path.exists(folder):
        os.makedirs(folder)

def create_single_tetrahedron(position, parity, r, deviation):
    # position: 3D array
    # parity: 1 (anti-clockwise) or -1 (clockwise)
    # r: 3D array
    # deviation: 6D array
    r1 = r[0]
    r2 = r[1]
    r3 = r[2]
    delta_theta1 = deviation[0]
    delta_theta2 = deviation[1]
    delta_theta3 = deviation[2]
    delta_r1 = deviation[3]
    delta_r2 = deviation[4]
    delta_r3 = deviation[5]
    
    primary = position

    transform = np.array([[1, 0, 0],
                          [0, np.cos(delta_theta1), np.sin(delta_theta1)],
                          [0, -np.sin(delta_theta1), np.cos(delta_theta1)]]) @ \
                        np.array([[np.cos(delta_theta2), 0, np.sin(delta_theta2)],
                                  [0, 1, 0],
                                  [-np.sin(delta_theta2), 0, np.cos(delta_theta2)]]) @ \
                        np.array([[np.cos(delta_theta3), np.sin(delta_theta3), 0],
                                  [-np.sin(delta_theta3), np.cos(delta_theta3), 0],
                                  [0, 0, 1]])

    r_initial = np.array([parity * (r1 + delta_r1), r2 + delta_r2, r3 + delta_r3])

    secondary1 = primary + transform @ np.array([r_initial[0], 0, 0])
    secondary2 = primary + transform @ np.array([0, r_initial[1], 0])
    secondary3 = primary + transform @ np.array([0, 0, r_initial[2]])
    return np.vstack((primary, secondary1, secondary2, secondary3))


def generate_random_deviations(deviation_range):
    deviation = np.array([random.uniform(-np.pi, np.pi) for i in range(3)] +
                         [random.uniform(deviation_range[i][0], deviation_range[i][1]) for i in range(3)])
    return deviation


def create_multiple_tetrahedra(vertices, parity, r, deviation_range):
    num_vert = vertices.shape[0]

    multiple_tetrahedra = []

    for i in range(num_vert):
        deviation = generate_random_deviations(deviation_range)
        single_tetrahedron = create_single_tetrahedron(vertices[i], parity, r, deviation)
        multiple_tetrahedra.append(single_tetrahedron)
    return np.array(multiple_tetrahedra)

def plot_tetrahedra(tetrahedra):
    global parity
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    for vertices in tetrahedra:
        # Plot the vertices as points
        ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], c="r", s=2)

        # Plot the edges of the tetrahedron
        edges = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
        for edge in edges:
            ax.plot([vertices[edge[0], 0], vertices[edge[1], 0]],
                    [vertices[edge[0], 1], vertices[edge[1], 1]],
                    [vertices[edge[0], 2], vertices[edge[1], 2]], "k-")

    # Set axis labels
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    # Set the aspect ratio to "equal" to make the tetrahedra appear in their true shape
    ax.set_box_aspect([1, 1, 1])

    # Show the plot
    plt.savefig(f"Figure/tetrahedra{parity}.png")
    plt.show()

# space = np.array([[0, 1000], [0, 1000], [0, 1000]])
# vertices = generate_3d_random(150, space, 60)
parity = 1
# r = [10, 20, 30]
# deviation_range = np.array([[0, 1], [-2, 2], [-1, 0]])
# tetrahedra = create_multiple_tetrahedra(vertices, parity, r, deviation_range)
# np.save(f"tetrahedra{parity}.npy", tetrahedra)
tetrahedra = np.load(f"tetrahedra{parity}.npy")
# plot_tetrahedra(tetrahedra)

radial_bins = []
for i in range(10):
    r1 = np.linspace(0, 1000, 10 + 1)[i]
    for j in range(i + 1, 10):
        r2 = np.linspace(0, 1000, 10 + 1)[j]
        if r2 - r1 >= 100:
            for k in range(j + 1, 10):
                r3 = np.linspace(0, 1000, 10 + 1)[k]
                if r3 - r2 >= 100:
                    radial_bins.append([r1, r2, r3])

def zeta_tetrahedra(l1, l2, l3, tetrahedra, radial_bins):
    global parity
    vertices = tetrahedra.reshape(-1, 3)
    zeta = []
    for bins in tqdm(radial_bins, desc="Processing bins"):
        bins_min = np.array(bins)
        bins_max = bins_min + int(1000 / 10)
        zeta_bin = np.imag(fe.estimator(l1, l2, l3, vertices, bins_min, bins_max, 1))
        zeta.append(zeta_bin)
    return zeta

zeta = zeta_tetrahedra(1, 1, 1, tetrahedra, radial_bins)
r_values = [r1 * r2 * r3 * z for (r1, r2, r3), z in zip(radial_bins, zeta)]
plt.plot(np.arange(len(radial_bins)), r_values)
plt.xlabel("Bin Index")
plt.ylabel(fr"$r_{1}r_{2}r_{3}\zeta(r_{1}r_{2}r_{3})$")
plt.title(r"$l_{1}=1$, $l_{2}=1$, $l_{3}=1$")
plt.savefig(f"Figure/zeta_tetrahedra{parity}.png")
plt.show()