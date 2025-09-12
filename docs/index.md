# popexposure

<div align="center">
  <img src="assets/popexposure_logo.png" alt="PopExposure Logo: The whole logo is a blue circle, and at the bottom it says 'popexposure' in white text. Above that, a yellow triangle is superimposed on a grid of nine white dots, representing a hazard and people being exposed to that hazard." width="200">
</div>

`popexposure` is an open-source Python package providing fast, memory-efficient, and consistent estimates of the number of people living near environmental hazards. It enables environmental scientists to assess population-level exposure to environmental hazards based on residential proximity.

!!! info "Methodological Details"
For comprehensive methodological details, see [McBrien et al (2025)]().

### Why `popexposure`?

Environmental epidemiologists often assess exposure to hazards using residential proximity (i.e., they consider an individual exposed if they live near a hazard). This computation presents technical difficulties, and different research teams usually apply their own solution. We developed popexposure to which quickly, efficiently, and consistently estimates the number of people living near environmental hazards.

- **Quick**: Optimized for processing large, fine-scale spatial datasets (e.g., exposure to oil and gas wells, which total millions of exposure points in the US) or datasets that cover a large area (e.g., national or global analyses of exposure).
- **Memory-efficient**: Only the necessary chunks of large population raster data are processed, avoiding loading the raster into memory.
- **Consistent**: `popexposure` implements a standardized methodology to ensure results are reproducible both within and across research teams.
- **Flexible**: `popexposure` can estimate the number of people exposesd to any type of hazard and according to any administrative boundary.

### Core API

| Function          | Overview                                                                                  | Inputs                                                                                      | Outputs                                                       |
| ----------------- | ----------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| `PopEstimator`    | Main class for estimating population exposure; initializes with population and admin data | `pop_data` (raster path), `admin_data` (GeoJSON/shapefile/GeoDataFrame, optional)           | PopEstimator object                                           |
| `est_exposed_pop` | Estimates number of people living within hazard buffer(s) using a population raster       | `hazard_data` (GeoJSON/shapefile/GeoDataFrame), `hazard_specific` (bool), `stat` (optional) | DataFrame with exposed population counts by hazard/admin unit |
| `est_total_pop`   | Estimates total population in administrative geographies using a population raster        | `stat` (optional, default `"sum"`)                                                          | DataFrame with total population per administrative unit       |

### What's Next?

- **[Installation](installation.md)** - Get PopExposure installed and ready to use
- **[Quick Start](quickstart.md)** - Learn the basics with a step-by-step tutorial
- **[Tutorials](tutorials/01_purpose_and_data_setup.ipynb)** - Explore real-world examples and advanced use cases
- **[API Reference](api/overview.md)** - Detailed documentation for all functions and classes

### Getting help and contributing

If you have any questions, a feature request, or would like to report a bug, please [open an issue](https://github.com/heathermcb/popexposure/issues). We also welcome any new contributions and ideas. If you want to add code, please submit a [pull request](https://github.com/heathermcb/popexposure/pulls) and we will get back to you when we can. Thanks!

### Citing this package

Please cite our paper [McBrien et al (2025)]().

### Authors

- [Heather McBrien](https://scholar.google.com/citations?user=0Hz3a1AAAAAJ&hl=en&oi=ao)
- [Joan A. Casey](https://scholar.google.com/citations?user=LjrwHBMAAAAJ&hl=en)
- [Lawrence Chillrud](https://scholar.google.com/citations?hl=en&user=HrSjGh0AAAAJ)
- [Nina M. Flores](https://scholar.google.com/citations?user=fkttN9UAAAAJ&hl=en&oi=ao)
- [Lauren B. Wilner](https://scholar.google.com/citations?user=rLX9LVYAAAAJ&hl=en&oi=ao)

### References

Our package is a fancy wrapper for the package [exactextract](https://pypi.org/project/exactextract/).
