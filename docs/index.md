# popexposure

`popexposure` is an open-source Python package providing fast, memory-efficient, and consistent estimates of the number of people living near environmental hazards. It enables environmental scientists to assess population-level exposure to environmental hazards based on residential proximity.

!!! info "Methodological Details"
    For comprehensive methodological details, see [McBrien et al (2025)]().

## Why `popexposure`?

Environmental epidemiologists often assess exposure to hazards using residential proximity (i.e., they consider an individual exposed if they live near a hazard). This computation presents technical difficulties, and different research teams usually apply their own solution. We developed an open-source Python package, popexposure, which quickly, efficiently, and consistently estimates the number of people living near environmental hazards. [LBW COMMENT: heather, i took this from the manuscript, feel free to edit]

- **Quick**: Optimized for processing large, fine-scale spatial datasets (e.g., exposure to oil and gas wells, which total millions of exposure points in the US) or datasets that cover a large area (e.g., national or global analyses of exposure)
- **Memory-efficient**: [LBW COMMENT: we never actually say how we did this in the paper -- heather, i will let you write in here what you want, and i think prob worth adding a sentence or two to the paper about this.]
- **Consistent**: `popexposure` implements a standardized methodology to ensure results are reproducible both within and across research teams.
- **Flexible**: `popexposure` can estimate the number of people exposesd to any type of hazard and according to any administrative boundary. 


## Available functions

| Function      | Overview                                                                 | Inputs                                                                                                      | Outputs                                                         |
| ------------- | ------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------- |
| `prep_data`   | Reads, cleans, and preprocesses geospatial hazard or admin unit data by removing empty or missing geometries, and buffering hazard data according to user-passed buffer distances    | Path to hazard or administrative unit file (`.geojson` or `.parquet`), `geo_type` (`"hazard"` or `"admin_unit"`) | Cleaned `GeoDataFrame` with valid geometries                    |
| `est_exposed_pop` | Estimates number of people living within hazard buffer(s) using a raster | Population raster path (`.tif`), hazard data, `hazard_specific` (bool), optional administrative units              | DataFrame with exposed population counts by hazard/administrative unit |
| `est_pop`         | Estimates total population in admin geographies using a raster           | Population raster path (`.tif`), administrative unit data (`GeoDataFrame`)                                         | DataFrame with total population per administrative unit                |

<!-- ## What's Next?

<div class="grid cards" markdown>

-   :material-download:{ .lg .middle } **Installation**

    ---

    Get PopExposure installed and ready to use

    [:octicons-arrow-right-24: Installation Guide](installation.md)

-   :material-rocket-launch:{ .lg .middle } **Quick Start**

    ---

    Learn the basics with a step-by-step tutorial

    [:octicons-arrow-right-24: Quick Start](quickstart.md)

-   :material-book-open:{ .lg .middle } **Tutorials**

    ---

    Explore real-world examples and advanced use cases

    [:octicons-arrow-right-24: View Tutorials](tutorials/basic_usage.md)

-   :material-api:{ .lg .middle } **API Reference**

    ---

    Detailed documentation for all functions and classes

    [:octicons-arrow-right-24: API Docs](api/overview.md)

</div> -->

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
