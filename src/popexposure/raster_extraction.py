"""
Raster extraction operations for population exposure analysis.

This module provides functionality for extracting values from raster datasets
using geometric masks, commonly used for population data extraction from
buffered hazard areas.
"""

import geopandas as gpd
import pandas as pd
import rasterio
from exactextract import exact_extract


class RasterExtractor:
    """Handles extraction of raster values for population exposure calculation."""

    @staticmethod
    def mask_raster_partial_pixel(
        shp_df: gpd.GeoDataFrame, raster_path: str, stat: str = "sum"
    ) -> gpd.GeoDataFrame:
        """
        Calculate the sum or mean of raster values (e.g., population) within each
        geometry using exact_extract for area-weighted extraction.

        Processes all geometry columns (main 'geometry' and any 'buffered_hazard_*'
        columns) and creates corresponding 'exposed' columns with raster value statistics.

        Can be used with a raster with any CRS.

        Parameters
        ----------
        shp_df : geopandas.GeoDataFrame
            Input GeoDataFrame containing geometries to extract raster values for.
            Must contain ID columns (ID_hazard, ID_admin_unit) and geometry columns.
        raster_path : str
            Path to the raster file (e.g., population raster)
        stat : str, default "sum"
            Statistic to calculate. Options: "sum" or "mean"

        Returns
        -------
        geopandas.GeoDataFrame
            Subset GeoDataFrame containing only ID columns and 'exposed' columns.
            Column naming: 'buffered_hazard_500' → 'exposed_500', 'geometry' → 'exposed'

        Notes
        -----
        - Can be used with a raster with any CRS
        - Uses exact_extract for area-weighted pixel extraction
        - Automatically reprojects geometries to match raster CRS
        - Invalid/empty geometries receive 0 values in exposed columns
        - Only returns ID and exposed columns, not original geometries
        - When stat="sum": returns total population within geometry
        - When stat="mean": returns average raster value within geometry

        Examples
        --------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>>
        >>> # Input data with buffered hazards
        >>> data = {
        ...     'ID_hazard': ['h1', 'h2'],
        ...     'geometry': [Point(0, 0).buffer(0.01), Point(1, 1).buffer(0.01)],
        ...     'buffered_hazard_500': [Point(0, 0).buffer(0.005), Point(1, 1).buffer(0.005)]
        ... }
        >>> gdf = gpd.GeoDataFrame(data, crs='EPSG:4326')
        >>>
        >>> # Extract total population (sum)
        >>> result_sum = RasterExtractor.mask_raster_partial_pixel(gdf, "population.tif", stat="sum")
        >>>
        >>> # Extract average population density (mean)
        >>> result_mean = RasterExtractor.mask_raster_partial_pixel(gdf, "population.tif", stat="mean")
        """
        with rasterio.open(raster_path) as src:
            raster_crs = src.crs

        geom_cols = [
            col
            for col in shp_df.columns
            if col.startswith("buffered_hazard") or col == "geometry"
        ]

        for geom_col in geom_cols:
            temp_gdf = shp_df[[geom_col]].copy()
            temp_gdf = temp_gdf.rename(columns={geom_col: "geometry"})
            temp_gdf = gpd.GeoDataFrame(temp_gdf, geometry="geometry", crs=shp_df.crs)

            if temp_gdf.crs != raster_crs:
                temp_gdf = temp_gdf.to_crs(raster_crs)

            # Identify invalid or empty geometries
            valid_mask = (
                temp_gdf.geometry.notnull()
                & temp_gdf.geometry.is_valid
                & (~temp_gdf.geometry.is_empty)
            )

            # Prepare a result column filled with zeros
            result = pd.Series(0, index=temp_gdf.index)

            # Only run exact_extract on valid geometries
            if valid_mask.any():
                valid_gdf = temp_gdf[valid_mask]
                num_exposed = exact_extract(
                    raster_path, valid_gdf, stat, output="pandas"
                )
                # Use appropriate data type based on statistic
                if stat == "sum":
                    result.loc[valid_mask] = num_exposed[stat].values.astype(int)
                else:  # stat == "mean"
                    result.loc[valid_mask] = num_exposed[stat].values.astype(float)

            exposed_col = f"exposed{geom_col.replace('buffered_hazard', '')}"
            shp_df[exposed_col] = result

        cols = [
            col
            for col in shp_df.columns
            if col.startswith("exposed") or col in ["ID_hazard", "ID_admin_unit"]
        ]
        shp_exposed = shp_df[cols]

        return shp_exposed
