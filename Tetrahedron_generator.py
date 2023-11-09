import numpy as np
import random
from random_number_generator import generate_1d_random, generate_3d_random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import fourPCF_estimator as fe

folders = ["Program/Figure"]
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
    plt.savefig(f"Program/Figure/tetrahedra_{parity}.png")
    plt.show()

space = np.array([[0, 1000], [0, 1000], [0, 1000]])
vertices = generate_3d_random(1500, space, 60)
parity = 1
r = [10, 20, 30]
deviation_range = np.array([[0, 1], [-2, 2], [-1, 0]])
tetrahedra = create_multiple_tetrahedra(vertices, parity, r, deviation_range)

plot_tetrahedra(tetrahedra)