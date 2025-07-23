# API Reference

PopExposure provides a simple, consistent API for population exposure analysis. The main entry point is the `PopEstimator` class.

## Core Classes and Functions

### PopEstimator Class

::: popexposure.find_exposure.PopEstimator

The `PopEstimator` class is the main interface for all population exposure calculations. It provides methods for data preparation and exposure estimation.

#### Key Methods

- **`prep_data()`**: Prepare and clean geospatial data
- **`est_exposed_pop()`**: Calculate population exposure to hazards  
- **`est_pop()`**: Estimate total population in administrative areas

### Function Reference

::: popexposure.find_exposure.PopEstimator.prep_data

::: popexposure.find_exposure.PopEstimator.exposed_pop

::: popexposure.find_exposure.PopEstimator.pop

## Data Requirements

### Input Data Formats

PopExposure accepts several common geospatial formats:

**Raster Data:**
- GeoTIFF (.tif, .tiff)
- Other GDAL-supported raster formats

**Vector Data:**
- GeoJSON (.geojson)
- Parquet with geometry (.parquet)
- Shapefile (.shp)
- Other OGR-supported vector formats

### Coordinate Reference Systems

- Input data should use the same coordinate reference system (CRS)
- Population rasters should use equal-area projections for accurate area calculations
- The package will attempt to handle CRS mismatches but results may be less accurate

## Error Handling

PopExposure includes comprehensive error checking for:

- Invalid file paths
- Mismatched coordinate systems
- Empty or invalid geometries
- Insufficient memory for large datasets

Common exceptions:

- `FileNotFoundError`: Input files cannot be found
- `ValueError`: Invalid parameters or data formats
- `MemoryError`: Insufficient memory for processing
