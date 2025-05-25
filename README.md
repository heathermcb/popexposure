<p align="left">
  <img src="writing/figs/modified-logo.png" alt="" width="120"/>
</p>

## popexposure: Functions to estimate the number of people living near environmental hazards

![Python](https://img.shields.io/badge/python-3.11-blue.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
[![GitHub](https://img.shields.io/badge/GitHub-Repo-black?logo=github)](https://github.com/heathermcb/popexposure)
[![PyPI version](https://badge.fury.io/py/popexposure.svg)](https://badge.fury.io/py/popexposure)

## Overview

`popexposure` is an open-source Python package providing fast, memory-efficient, and consistent estimates of the number of people living near environmental hazards, enabling environmental scientists to assess population-level exposure to environmental hazards based on residential proximity. Methodological details can be found in [McBrien et al (2025)](). Extensive documentation can be found on in our quick start [tutorial](https://github.com/heathermcb/popexposure/blob/main/demo/).

## Installation

The easiest way to install `popexposure` is via the latest pre-compiled binaries from PyPI with:

```bash
pip install popexposure
```

You can build `popexposure` from source as you would any other Python package with:

```bash
git clone https://github.com/heathermcb/popexposure
cd popexposure
python -m pip install .
```

## Tutorials

A number of tutorials providing worked examples using `popexposure` can be found in our [demos](https://github.com/heathermcb/Pop_Exp/tree/main/demo/demo) folder.

## Quickstart

```python
import glob
import pandas as pd
from popexposure.find_exposure import PopEstimator

# Instantiate estimator
pop_est = PopEstimator()

# List of years and corresponding hazard file paths
years = [2016, 2017, 2018]
hazard_paths = [
    "hazard_2016.geojson",
    "hazard_2017.geojson",
    "hazard_2018.geojson"
]
pop_path = "my_pop_raster.tif"
admin_units_path = "my_admin_units.geojson"

# Prepare admin units data
admin_units = pop_est.prep_data(admin_units_path, geo_type="admin_unit")

# Find total num ppl residing <= 10km of each hazard in each year
exposed_list = []

for year, hazard_path in zip(years, hazard_paths):
    # Prepare hazard data for this year
    hazards = pop_est.prep_data(hazard_path, geo_type="hazard")
    # Estimate exposed population
    exposed = pop_est.est_exposed_pop(
        pop_path=pop_path,
        hazard_specific=False,  # set to True if you want per-hazard results
        hazards=hazards,
        admin_units=admin_units # leave as None for estimates not by admin unit
    )
    exposed['year'] = year
    exposed_list.append(exposed)

exposed_df = pd.concat(exposed_list, axis=0)

# Save output
exposed_df.to_parquet("pop_exposed_to_hazards.parquet")
```

## Available functions

| Function      | Overview                                                                 | Inputs                                                                                                      | Outputs                                                         |
| ------------- | ------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------- |
| `prep_data`   | Reads, cleans, and preprocesses geospatial hazard or admin unit data by removing empty or missing geometries, and buffering hazard data according to user-passed buffer distances    | Path to hazard or administrative unit file (`.geojson` or `.parquet`), `geo_type` (`"hazard"` or `"admin_unit"`) | Cleaned `GeoDataFrame` with valid geometries                    |
| `est_exposed_pop` | Estimates number of people living within hazard buffer(s) using a raster | Population raster path (`.tif`), hazard data, `hazard_specific` (bool), optional administrative units              | DataFrame with exposed population counts by hazard/administrative unit |
| `est_pop`         | Estimates total population in admin geographies using a raster           | Population raster path (`.tif`), administrative unit data (`GeoDataFrame`)                                         | DataFrame with total population per administrative unit                |

## Getting help and contributing

If you have any questions, a feature request, or would like to report a bug, please [open an issue](https://github.com/heathermcb/Pop_Exp/issues). We also welcome any new contributions and ideas. If you want to add code, please submit a [pull request](https://github.com/heathermcb/Pop_Exp/pulls) and we will get back to you when we can. Thanks!

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
