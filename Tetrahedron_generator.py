import numpy as np
import random

def create_single_tetrahedron(position, parity, r, deviation):
    transform_matrix = np.array([[1, 0, 0],
                                 [0, np.cos(deviation[0]), np.sin(deviation[0])],
                                 [0, -np.sin(deviation[0]), np.cos(deviation[0])]]) @ \
                       np.array([[np.cos(deviation[1]), 0, np.sin(deviation[1])],
                                 [0, 1, 0],
                                 [-np.sin(deviation[1]), 0, np.cos(deviation[1])]]) @ \
                       np.array([[np.cos(deviation[2]), np.sin(deviation[2]), 0],
                                 [-np.sin(deviation[2]), np.cos(deviation[2]), 0],
                                 [0, 0, 1]])

    r_before_trans = np.array([parity * (r[0] + deviation[3]), r[1] + deviation[4], r[2] + deviation[5]])

    primary = position
    secondary_1 = primary + transform_matrix @ np.array([r_before_trans[0], 0, 0])
    secondary_2 = primary + transform_matrix @ np.array([0, r_before_trans[1], 0])
    secondary_3 = primary + transform_matrix @ np.array([0, 0, r_before_trans[2]])

    return np.vstack((primary, secondary_1, secondary_2, secondary_3))


def generate_random_deviations(r_min, r_max):
    return np.array([random.uniform(-np.pi, np.pi) for _ in range(3)] +
                    [random.uniform(r_min[i], r_max[i]) for i in range(3)])


def create_multiple_tetrahedra(vertices, parity, r, deviation_range):
    num_of_vertices = vertices.shape[0]
    deviation_min, deviation_max = deviation_range[:3], deviation_range[3:]

    multiple_tetrahedra = []

    for i in range(num_of_vertices):
        deviation = generate_random_deviations(deviation_min, deviation_max)
        single_tetrahedron = create_single_tetrahedron(vertices[i], parity, r, deviation)
        multiple_tetrahedra.append(single_tetrahedron)

    return np.array(multiple_tetrahedra)