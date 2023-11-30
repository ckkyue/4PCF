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
# vertices = generate_3d_random(1500, space, 60)
# parity = 1 # 1 or -1
# r = [10, 20, 30]
# deviation_range = np.array([[0, 1], [-2, 2], [-1, 0]])
# tetrahedra = create_multiple_tetrahedra(vertices, parity, r, deviation_range)
# np.save(f"tetrahedra{parity}.npy", tetrahedra)
# plot_tetrahedra(tetrahedra)

def zeta_tetrahedra(l1, l2, l3, tetrahedra, bins_min, bins_max):
    global parity
    weights = [1 / len(tetrahedra)]*len(tetrahedra)
    zeta = []
    for tetrahedron in tqdm(tetrahedra):
        vertices = tetrahedron
        zeta_tetrahedron = np.imag(fe.estimator(l1, l2, l3, vertices, bins_min, bins_max, weights))
        zeta.append(zeta_tetrahedron)
    return zeta

tetrahedra_ccw = np.load("tetrahedra1.npy")
tetrahedra_cw = np.load("tetrahedra-1.npy")

# zeta111_ccw = zeta_tetrahedra(1, 1, 1, tetrahedra_ccw, [5, 15, 25], [15, 25, 35])
# np.save("zeta111_tetrahedra1.npy", zeta111_ccw)
# zeta111_cw = zeta_tetrahedra(1, 1, 1, tetrahedra_cw, [5, 15, 25], [15, 25, 35])
# np.save("zeta111_tetrahedra-1.npy", zeta111_cw)
zeta111_ccw = np.load("zeta111_tetrahedra1.npy")
zeta111_cw = np.load("zeta111_tetrahedra-1.npy")
plt.plot(np.arange(len(zeta111_ccw)), zeta111_ccw, color="blue", label=r"$r_{1}\in[5, 15]$, $r_{2}\in[15, 25]$, $r_{3}\in[25, 35]$ (ccw)")
plt.plot(np.arange(len(zeta111_cw)), zeta111_cw, color="orange", label=r"$r_{1}\in[5, 15]$, $r_{2}\in[15, 25]$, $r_{3}\in[25, 35]$ (cw)")
plt.axhline(np.mean(zeta111_ccw), linestyle="--", color="black", alpha=0.5)
plt.axhline(np.mean(zeta111_cw), linestyle="--", color="black", alpha=0.5)
plt.xlabel("Tetrahedron Index")
plt.ylabel(r"$\zeta_{l_{1}, l_{2}, l_{3}}(r_{1}, r_{2}, r_{3})$")
plt.title(r"$l_{1}=1$, $l_{2}=1$, $l_{3}=1$")
plt.legend()
plt.savefig(f"Figure/zeta111_tetrahedra_mixed.png")
plt.show()

# zeta122_ccw = zeta_tetrahedra(1, 2, 2, tetrahedra_ccw, [5, 15, 25], [15, 25, 35])
# np.save("zeta122_tetrahedra1.npy", zeta122_ccw)
# zeta122_cw = zeta_tetrahedra(1, 2, 2, tetrahedra_cw, [5, 15, 25], [15, 25, 35])
# np.save("zeta122_tetrahedra-1.npy", zeta122_cw)
zeta122_ccw = np.load("zeta122_tetrahedra1.npy")
zeta122_cw = np.load("zeta122_tetrahedra-1.npy")
plt.plot(np.arange(len(zeta122_ccw)), zeta122_ccw, color="blue", label=r"$r_{1}\in[5, 15]$, $r_{2}\in[15, 25]$, $r_{3}\in[25, 35]$ (ccw)")
plt.plot(np.arange(len(zeta122_cw)), zeta122_cw, color="orange", label=r"$r_{1}\in[5, 15]$, $r_{2}\in[15, 25]$, $r_{3}\in[25, 35]$ (cw)")
plt.axhline(np.mean(zeta122_ccw), linestyle="--", color="black", alpha=0.5)
plt.axhline(np.mean(zeta122_cw), linestyle="--", color="black", alpha=0.5)
plt.xlabel("Tetrahedron Index")
plt.ylabel(r"$\zeta_{l_{1}, l_{2}, l_{3}}(r_{1}, r_{2}, r_{3})$")
plt.title(r"$l_{1}=1$, $l_{2}=2$, $l_{3}=2$")
plt.legend()
plt.savefig(f"Figure/zeta122_tetrahedra_mixed.png")
plt.show()