# popexposure <a href="https://heathermcb.github.io/popexposure/"><img src="docs/assets/popexposure_logo.png" align="right" alt="popexposure documentation website" width="120" /></a>

![Python](https://img.shields.io/badge/python-3.11-blue.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
[![GitHub](https://img.shields.io/badge/GitHub-Repo-black?logo=github)](https://github.com/heathermcb/popexposure)
![Coverage](https://img.shields.io/badge/coverage-100%25-brightgreen)
[![PyPI version](https://badge.fury.io/py/popexposure.svg)](https://badge.fury.io/py/popexposure)

## Overview

`popexposure` is an open-source Python package providing fast, memory-efficient, and consistent estimates of the number of people living near environmental hazards, enabling environmental scientists to assess population-level exposure to environmental hazards based on residential proximity. Methodological details can be found in [McBrien et al (2025)](). Extensive documentation can be found on our [website](https://heathermcb.github.io/popexposure/), interactive examples in our [tutorials](https://github.com/heathermcb/popexposure/tree/main/docs/tutorials), or scaffolding code in our [quickstart](https://heathermcb.github.io/popexposure/quickstart/).

## Installation

The easiest way to install `popexposure` is via the latest pre-compiled binaries from PyPI with:

```bash
pip install popexposure
```

You can build `popexposure` from source with:

```bash
git clone https://github.com/heathermcb/popexposure
cd popexposure
python -m pip install .
```

## Documentation and tutorials

The `popexposure` API comes with a comprehensive documentation [website](https://heathermcb.github.io/popexposure/).

The [tutorials](https://github.com/heathermcb/popexposure/tree/main/docs/tutorials) folder contains three jupyter notebooks (and a data download bash script) than can be run in sequence to interactively learn how to use `popexposure` to estimate the number of people residing near California wildfire disasters for each year spanning 2016 to 2018.

## Quickstart

```python
import popexposure as ex
import geopandas as gpd
import pandas as pd

# Set paths
my_pop_raster_path = "path/to/my_pop_raster.tif"
my_admin_units_path = "path/to/my_admin_units.geojson"
my_hazard_data_path = "path/to/my_hazards.geojson"

# Read in hazard data and define buffer(s) of interest (in meters)
hazard_gdf = gpd.read_file(my_hazard_data_path)
hazard_gdf["buffer_dist_0km"] = 0
hazard_gdf["buffer_dist_10km"] = 10_000
hazard_gdf["buffer_dist_20km"] = 20_000

# Instantiate estimator
pop_est = ex.PopEstimator(pop_data=my_pop_raster_path, admin_data=my_admin_units_path)

# Find total num ppl residing <= 0, 10, 20km of each hazard in the hazard_gdf
exposed_df = pop_est.est_exposed_pop(
    hazard_data=hazard_gdf,  # can also pass in a filepath here, assuming it has necessary columns
    hazard_specific=False,   # set to True if you need per-hazard results, False for cumulative exposure
)

# Save output
exposed_df.to_parquet("pop_exposed_to_hazards.parquet")
```

## Available methods

| Function          | Overview                                                                                           | Inputs                                                             | Outputs                                                       |
| ----------------- | -------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ | ------------------------------------------------------------- |
| `PopEstimator`    | Main class for estimating population exposure; initializes with population and optional admin data | `pop_data` (raster path), `admin_data` (GeoJSON or shapefile path) | PopEstimator object                                           |
| `est_exposed_pop` | Estimates number of people living within a specified distance of hazards                           | `hazards` (GeoJSON/shapefile), `hazard_specific` (bool)            | DataFrame with exposed population counts by hazard/admin unit |
| `est_total_pop`   | Estimates total population in administrative units                                                 | None (uses data provided at initialization)                        | DataFrame with total population per administrative unit       |

## Getting help and contributing

If you have any questions, a feature request, or would like to report a bug, please [open an issue](https://github.com/heathermcb/popexposure/issues). We also welcome any new contributions and ideas. If you want to add code, please submit a [pull request](https://github.com/heathermcb/popexposure/pulls) and we will get back to you when we can. Thanks!

## Citing this package

Please cite our paper [McBrien et al (2025)]().

## Authors

- [Heather McBrien](https://scholar.google.com/citations?user=0Hz3a1AAAAAJ&hl=en&oi=ao)
- [Joan A. Casey](https://scholar.google.com/citations?user=LjrwHBMAAAAJ&hl=en)
- [Lawrence Chillrud](https://scholar.google.com/citations?hl=en&user=HrSjGh0AAAAJ)
- [Nina M. Flores](https://scholar.google.com/citations?user=fkttN9UAAAAJ&hl=en&oi=ao)
- [Lauren B. Wilner](https://scholar.google.com/citations?user=rLX9LVYAAAAJ&hl=en&oi=ao)

## References

Our package is a fancy wrapper for the package [exactextract](https://pypi.org/project/exactextract/).
