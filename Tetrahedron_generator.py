import numpy as np
import random

def create_single_tetrahedron(position, parity, r, deviation):
    # position: 3D array

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


def generate_random_deviations(r_min, r_max):
    return np.array([random.uniform(-np.pi, np.pi) for i in range(3)] +
                    [random.uniform(r_min[i], r_max[i]) for i in range(3)])


def create_multiple_tetrahedra(vertices, parity, r, deviation_range):
    num_vert = vertices.shape[0]
    deviation_min, deviation_max = deviation_range[:3], deviation_range[3:]

    multiple_tetrahedra = []

    for i in range(num_vert):
        deviation = generate_random_deviations(deviation_min, deviation_max)
        single_tetrahedron = create_single_tetrahedron(vertices[i], parity, r, deviation)
        multiple_tetrahedra.append(single_tetrahedron)

    return np.array(multiple_tetrahedra)