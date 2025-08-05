# API Reference

PopExposure provides a simple, consistent API for population exposure analysis. The main entry point is the `PopEstimator` class.

## Core Classes and Functions

### PopEstimator Class

::: popexposure.find_exposure.PopEstimator

The `PopEstimator` class is the main interface for all population exposure calculations. It provides methods for data preparation and exposure estimation.

#### Key Methods

- **`est_exposed_pop()`**: Calculate population exposure to hazards
- **`est_total_pop()`**: Estimate total population in administrative areas

### Function Reference

::: popexposure.find_exposure.PopEstimator.est_exposed_pop

::: popexposure.find_exposure.PopEstimator.est_total_pop

## Data Requirements

### Input Data Formats

**Data Requirements:**

    - **Hazard data**: GeoJSON or GeoParquet, must contain string ``ID_hazard`` column with unique hazard IDs, ``buffer_dist_*`` numeric columns, and ``geometry`` column with geometry objects.
    - **Admin units**: GeoJSON or GeoParquet, must contain string ``ID_admin_unit``  column with unique admin IDs, and ``geometry`` column column with geometry objects.
    - **Population raster**: Any CRS supported

See tutorials and API docs for more detailed information on data requirements.
