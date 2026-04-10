# 4PCF

This repository contains a small research codebase for experimenting with the four-point correlation function (4PCF) and parity-odd modes in galaxy catalogs and toy tetrahedron samples.

The current workflow is centered on three pieces:

- computing a 4PCF estimator from Cartesian point catalogs,
- generating toy tetrahedron catalogues to test parity sensitivity,
- converting BOSS-style redshift data into Cartesian positions and weights.

`edge_correction.py` is still experimental and is not yet integrated into the main estimator pipeline.

If you want the scientific background, see the final report in the repository and the references listed below.

## Repository Layout

### Core scripts

- [fourPCF_estimator.py](fourPCF_estimator.py) computes the 4PCF estimator for a set of vertices and radial bins.
- [tetrahedron_generator.py](tetrahedron_generator.py) builds and evaluates a toy tetrahedron model.
- [random_number_generator.py](random_number_generator.py) creates 1D and 3D random samples with a minimum separation.
- [read_data.py](read_data.py) reads raw catalog data and converts it into weights and Cartesian coordinates.
- [redshift_distance.py](redshift_distance.py) converts redshift into a comoving distance estimate.
- [edge_correction.py](edge_correction.py) contains work-in-progress coupling and edge-correction utilities.

### Sample data and outputs

- [test_sample.txt](test_sample.txt) is a small BOSS CMASS-style sample used for testing.
- [vertices_sample.npy](vertices_sample.npy) and [weights_sample.npy](weights_sample.npy) are precomputed Cartesian positions and weights for that sample.
- [tetrahedra1.npy](tetrahedra1.npy) and [tetrahedra-1.npy](tetrahedra-1.npy) store toy tetrahedra with positive and negative parity.
- [zeta111_tetrahedra1.npy](zeta111_tetrahedra1.npy), [zeta111_tetrahedra-1.npy](zeta111_tetrahedra-1.npy), [zeta122_tetrahedra1.npy](zeta122_tetrahedra1.npy), and [zeta122_tetrahedra-1.npy](zeta122_tetrahedra-1.npy) are saved estimator outputs.
- [Figure/](Figure/) is used for generated plots.

## Requirements

Install the Python dependencies with:

```bash
pip install numpy scipy sympy matplotlib tqdm
```

Using a virtual environment is recommended if you want to keep the scientific stack isolated from your system Python.

## Quick Start

The repository is set up around runnable scripts rather than a packaged Python module.

1. Run the main estimator example:

```bash
python fourPCF_estimator.py
```

This script loads [vertices_sample.npy](vertices_sample.npy) and [weights_sample.npy](weights_sample.npy), builds a radial-bin grid, and evaluates the 4PCF estimator for the selected multipole tuple.

2. Run the toy tetrahedron example:

```bash
python tetrahedron_generator.py
```

This script loads the saved tetrahedron samples, evaluates the toy-model 4PCF, and writes plots into [Figure/](Figure/).

3. If you want to regenerate the sample arrays from raw catalog data, use [read_data.py](read_data.py) together with [test_sample.txt](test_sample.txt).

## Main Functions

### [fourPCF_estimator.py](fourPCF_estimator.py)

- `longitude(x, y)` returns angles in the range $[0, 2\pi)$.
- `cart_to_sphe(cart)` converts Cartesian coordinates to spherical coordinates.
- `shell_vol(bin_min, bin_max)` computes the volume of a spherical shell.
- `a(l, m, primary_vert, secondary_vert, bin_min, bin_max, weights)` evaluates the spherical-harmonic coefficient used by the estimator.
- `estimator(l1, l2, l3, vertices, bins_min, bins_max, weights)` returns the 4PCF estimator for the selected multipoles and radial bins.

### [tetrahedron_generator.py](tetrahedron_generator.py)

- `create_single_tetrahedron(position, parity, r, deviation)` creates one tetrahedron from a seed position.
- `generate_random_deviations(deviation_range)` samples random angular and radial perturbations.
- `create_multiple_tetrahedra(vertices, parity, r, deviation_range)` expands a vertex catalog into a tetrahedron ensemble.
- `plot_tetrahedra(tetrahedra)` saves a 3D plot of the toy tetrahedra.
- `zeta_tetrahedra(l1, l2, l3, tetrahedra, bins_min, bins_max)` evaluates the 4PCF estimator on the toy ensemble.

### [read_data.py](read_data.py)

- `read(data_path)` loads the comma-separated input file.
- `get_weights(data)` computes the combined FKP/systematic/completeness weight used by the catalog.
- `get_carts(data)` converts right ascension, declination, and redshift into Cartesian coordinates.

### [redshift_distance.py](redshift_distance.py)

- `redshift_to_dist(z, type="DCMR", h=H0/100.0, Omega_m=Omega_m, n=1000)` converts redshift to a comoving distance estimate.

## Notes

- The estimator scripts currently read the `.npy` sample files directly at import time, so they are easiest to use as standalone scripts.
- `edge_correction.py` references estimator logic that is still under development, so treat it as experimental.
- The repository currently includes precomputed outputs so you can inspect results without rerunning the full pipeline.

## References

This project is motivated by parity-odd 4PCF measurements in BOSS galaxies, especially:

- Philcox (2022): https://journals.aps.org/prd/abstract/10.1103/PhysRevD.106.063501
- Hou et al. (2023): https://academic.oup.com/mnras/article/522/4/5701/7169316

BibTeX entries:

```bibtex
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
    year = {2023},
    month = {may},
    publisher = {Oxford University Press},
    volume = {522},
    number = {4},
    pages = {5701--5739},
    author = {Jiamin Hou and Zachary Slepian and Robert N Cahn},
    title = {{Measurement of parity-odd modes in the large-scale 4-point correlation function of Sloan Digital Sky Survey Baryon Oscillation Spectroscopic Survey twelfth data release CMASS and LOWZ galaxies}},
    journal = {Monthly Notices of the Royal Astronomical Society}
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
