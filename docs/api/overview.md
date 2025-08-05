# API Reference

**popexposure** provides simple, consistent methods for population exposure analysis. Users can access this through the `PopEstimator` class, which provides methods that prepare data and perform exposure estimation.

**Key Methods**

- **`est_exposed_pop()`**: Calculate population exposure to hazards.
- **`est_total_pop()`**: Estimate total population in administrative areas.

**Data Requirements:**

- **Hazard data**: `GeoJSON` or `GeoParquet`, must contain string `ID_hazard` column with unique hazard IDs, `buffer_dist_*` numeric columns, and `geometry` column with geometry objects.

- **Admin units**: `GeoJSON` or `GeoParquet`, must contain string `ID_admin_unit` column with unique admin IDs, and `geometry` column with geometry objects.

- **Population raster**: Any format supported by rasterio with any CRS.

See tutorials and API docs below for more detailed information on data requirements and examples.

::: popexposure.pop_estimator.PopEstimator
