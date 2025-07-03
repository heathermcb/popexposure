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
import functools


class RasterExtractor:
    """Handles extraction of raster values for population exposure calculation."""

    @staticmethod
    def mask_raster_partial_pixel(
        shp_df: gpd.GeoDataFrame, raster_path: str
    ) -> gpd.GeoDataFrame:
        """
        Calculate the sum of raster values (e.g., population) within each
        geometry in buffered hazard columns using exact_extract.

        Parameters
        ----------
        shp_df : geopandas.GeoDataFrame
            Input GeoDataFrame containing geometries to extract raster values for
        raster_path : str
            Path to the raster file

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame with 'exposed' columns containing the sum for each geometry
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
                    raster_path, valid_gdf, "sum", output="pandas"
                )
                result.loc[valid_mask] = num_exposed["sum"].values

            exposed_col = f"exposed{geom_col.replace('buffered_hazard', '')}"
            shp_df[exposed_col] = result

        cols = [
            col
            for col in shp_df.columns
            if col.startswith("exposed") or col in ["ID_hazard", "ID_admin_unit"]
        ]
        shp_exposed = shp_df[cols]

        return shp_exposed
