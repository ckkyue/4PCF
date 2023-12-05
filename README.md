# 4PCF
This repository contains the code for my recent project on cosmological parity violation. The primary code file, "fourPCF_estimator.py," is responsible for computing the 4PCF estimator for a set of discrete Cartesian points.
The "edge_correction.py" is currently under development and has not yet been incorporated into the calculation of the 4PCF estimator.

## fourPCF_estimator.py
This file contains the main code for calculating the 4PCF estimator.

### def longitude(x, y)
This function returns the longitude of a point within the range $[0, 2\pi)$, given two 1D arrays.

### def cart_to_sphe(cart)
This function converts Cartesian coordinates to spherical coordinates.

### def bin_func(x, bin_min, bin_max)
This function returns $1$ if $x$ falls within a specified bin range and $0$ otherwise.

### def shell_vol(bin_min, bin_max)
This function returns the volume of a spherical shell given two radii.

### def a(l, m, primary_vert, secondary_vert, bin_min, bin_max, weight)
This function is the core function for calculating the $a$ constant.

### def estimator(l1, l2, l3, vertices, bins_min, bins_max, weights)
This function returns the 4PCF estimator, given the angular momenta multipole `$(l_{1}, l_{2}, l_{3})$`.

## random_number_generator.py
This file contains the code for generating random points.

### min_sep_test(point, points, min_sep)
This function checks if the `point` is separated by at least `min_sep` distance from all points in the `points` collection.

### generate_1d_random(num, space, min_sep)
This function generates random 1D points with a minimum separation distance `min_sep`.

### generate_3d_random(num, space, min_sep)
This function generates random 3D points with a minimum separation distance `min_sep`.

## tetrahedron_generator.py
This file contains the code for generating the toy model of tetrahedra.

### def create_single_tetrahedron(position, parity, r, deviation)
This function returns the vertices of a single tetrahedron randomly.

### def generate_random_deviations(deviation_range)
This function generates random deviations in rotations of vertices and random distances, which will be added to the three sides of a tetrahedron.

### def create_multiple_tetrahedra(vertices, parity, r, deviation_range)
This function returns the vertices of multiple tetrahedra randomly, given the primary vertices.

### def plot_tetrahedra(tetrahedra)
This function plots the tetrahedra.

### def zeta_tetrahedra(l1, l2, l3, tetrahedra, bins_min, bins_max)
This function returns the 4PCF estimator of the toy model of tetrahedra, given the angular momenta multipole `$(l_{1}, l_{2}, l_{3})$`.

## tetrahedra1.npy
This file contains 1500 sample counterclockwise tetrahedra with a parity of $1$.

## tetrahedra-1.npy
This file contains 1500 sample clockwise tetrahedra with a parity of $-1$.

## test_sample.txt
This file contains the raw data of 1000 galaxies from the CMASS sample of BOSS.

## read_data.py
This file contains functions for reading and processing data.

### def read(data_path)
This function returns the data given the path to the raw data.

### def get_weights(data)
This function returns the weights given the data.

### def get_carts(data)
This function returns the Cartesian coordinates given the data.

## redshift_distance.py
This file contains the code to calculate the physical distance given the redshfit $z$.

### def redshift_to_dist(z, type="DCMR", h=H0/100.0, Omega_m=Omega_m, n=1000)
This function returns the physical distance given the redshift $z$.

## vertices_sample.npy
This file contains the Cartesian coordinates of 1000 galaxies from the CMASS sample of BOSS.

## weights_sample.npy
This file contains the weights associated with the 1000 galaxies from the CMASS sample of BOSS.
