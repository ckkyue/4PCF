# 4PCF
This repository contains the code for my recent project on cosmological parity violation. The primary code file, "fourPCF_estimator.py," is responsible for computing the 4PCF estimator for a set of discrete Cartesian points.
The "edge_correction.py" is currently under development and has not yet been incorporated into the calculation of the 4PCF estimator.

## fourPCF_estimator.py
This file contains the main code for calculating the 4PCF estimator.

### def longitude(x, y)
This function returns the longitude of a point within the range $[0, 2\pi)$, given two 1D arrays `x` and `y`.

### def cart_to_sphe(cart)
This function converts Cartesian coordinates `cart` to spherical coordinates.

### def shell_vol(bin_min, bin_max)
This function returns the volume of a spherical shell given two radii `bin_min` and `bin_max`.

### def a(l, m, primary_vert, secondary_vert, bin_min, bin_max, weight)
This function is the core function for calculating the $a$ coefficient.

### def estimator(l1, l2, l3, vertices, bins_min, bins_max, weights)
This function returns the 4PCF estimator in each bin `(bin_min, bin_max)`, given the angular momentum multipole $(l_{1}, l_{2}, l_{3})$.

## random_number_generator.py
This file contains the code for generating random points.

### min_sep_test(point, points, min_sep)
This function checks if the `point` is separated by at least `min_sep` distance from all points in the `points` collection.

### generate_1d_random(num, space, min_sep)
This function generates `num` random 1D points with a minimum separation distance `min_sep`.

### generate_3d_random(num, space, min_sep)
This function generates `num` random 3D points with a minimum separation distance `min_sep`.

## tetrahedron_generator.py
This file contains the code for generating the toy model of tetrahedra.

### def create_single_tetrahedron(position, parity, r, deviation)
This function returns the vertices of a single tetrahedron with random deviations.

### def generate_random_deviations(deviation_range)
This function generates random deviations in rotations of vertices and random distances, which will be added to the three sides of a tetrahedron.

### def create_multiple_tetrahedra(vertices, parity, r, deviation_range)
This function returns the vertices of multiple tetrahedra with random deviations, given the primary vertices.

### def plot_tetrahedra(tetrahedra)
This function creates a 3D plot of the `tetrahedra`.

### def zeta_tetrahedra(l1, l2, l3, tetrahedra, bins_min, bins_max)
This function returns the 4PCF estimator of the toy model of tetrahedra in each bin `(bin_min, bin_max)`, given the angular momentum multipole $(l_{1}, l_{2}, l_{3})$.

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
This function returns the weights given the `data`.

### def get_carts(data)
This function returns the Cartesian coordinates given the `data`.

## redshift_distance.py
This file contains the code to calculate the physical distance given the redshfit `z`.

### def redshift_to_dist(z, type="DCMR", h=H0/100.0, Omega_m=Omega_m, n=1000)
This function returns the physical distance given the redshift `z`.

## vertices_sample.npy
This file contains the Cartesian coordinates of 1000 galaxies from the CMASS sample of BOSS.

## weights_sample.npy
This file contains the weights associated with the 1000 galaxies from the CMASS sample of BOSS.

## Libraries
Ensure that you have installed the following libraries.
```
pip install matplotlib
pip install numpy
pip install scipy
pip install sympy
pip install tqdm
```

## References
The primary objective of this work is to verify the detection of parity-odd 4PCF by Philcox (2022) and Hou et al. (2023).

Philcox (2022): https://journals.aps.org/prd/abstract/10.1103/PhysRevD.106.063501

Hou et al. (2023): https://academic.oup.com/mnras/article/522/4/5701/7169316?login=true

To cite their papers, use the following BibTex:
```
@article{PhysRevD.106.063501,
    title = {{Probing parity violation with the four-point correlation function of BOSS galaxies}},
    author = {Philcox, Oliver H. E.},
    journal = {Phys. Rev. D},
    volume = {106},
    issue = {6},
    pages = {063501},
    numpages = {29},
    year = {2022},
    month = {Sep},
    publisher = {American Physical Society},
    doi = {10.1103/PhysRevD.106.063501},
    url = {https://link.aps.org/doi/10.1103/PhysRevD.106.063501}
}
@article{Hou_2023,
    doi = {10.1093/mnras/stad1062},
    url = {https://doi.org/10.1093%2Fmnras%2Fstad1062},
    year = 2023,
    month = {may},
    publisher = {Oxford University Press ({OUP})},
    volume = {522},
    number = {4},
    pages = {5701--5739},
    author = {Jiamin Hou and Zachary Slepian and Robert N Cahn},
    title = {{Measurement of parity-odd modes in the large-scale 4-point correlation function of Sloan Digital Sky Survey Baryon Oscillation Spectroscopic Survey twelfth data release {CMASS} and {LOWZ} galaxies}},
    journal = {Monthly Notices of the Royal Astronomical Society}
}
```
