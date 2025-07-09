"""
Spatial geometry operations for population exposure analysis.

This module provides classes for buffering, unioning, and intersection
operations used in population exposure calculations.
"""

import geopandas as gpd
import pandas as pd
import pyproj
from shapely.ops import unary_union, transform
from shapely.geometry.base import BaseGeometry
import functools
from .geometry_validator import GeometryValidator


class GeometryOperations:
    """Handles buffering, unioning, and intersection operations."""

    @staticmethod
    def get_buffered_geometry(row: pd.Series, buffer_col: str) -> BaseGeometry:
        """
        Buffer a single geometry using its appropriate UTM projection.

        Projects the geometry from WGS84 to UTM for accurate distance-based
        buffering, applies the buffer, then reprojects back to WGS84.

        Parameters
        ----------
        row : pandas.Series
            Row from a GeoDataFrame containing 'utm_projection', 'geometry',
            and the specified buffer distance column
        buffer_col : str
            Name of the column containing the buffer distance in meters

        Returns
        -------
        shapely.geometry
            Buffered geometry in WGS84 coordinate system

        Notes
        -----
        - Buffer distance is in meters
        - Uses UTM projection for accurate metric buffering
        - Transformers are created fresh for each call to ensure thread safety
        """
        # Set up transformers
        best_utm = row["utm_projection"]
        hazard_geom = row["geometry"]
        buffer_dist = row[buffer_col]

        # Set up transformers only once per call
        to_utm = pyproj.Transformer.from_crs(
            "EPSG:4326", best_utm, always_xy=True
        ).transform
        to_wgs = pyproj.Transformer.from_crs(
            best_utm, "EPSG:4326", always_xy=True
        ).transform

        # Project to UTM, buffer, then project back
        geom_utm = transform(to_utm, hazard_geom)
        buffered = geom_utm.buffer(buffer_dist)
        buffered_wgs = transform(to_wgs, buffered)
        return buffered_wgs

    @staticmethod
    def add_buffered_geometry_columns(
        gdf: gpd.GeoDataFrame,
    ) -> gpd.GeoDataFrame:
        """
        Add buffered geometry columns based on buffer_dist columns.

        For each column starting with 'buffer_dist_', creates a corresponding
        'buffered_hazard_' column with geometries buffered by the specified distance.
        Requires that the GeoDataFrame already has 'utm_projection' column.

        Each geometry is projected to its optimal UTM coordinate system for accurate
        distance-based buffering, then reprojected back to the original CRS.

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Input GeoDataFrame with buffer_dist columns and utm_projection column

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame with added buffered geometry columns

        Notes
        -----
        - Buffer distances are in meters
        - Column naming: 'buffer_dist_500' → 'buffered_hazard_500'
        - Uses UTM projection for accurate metric buffering
        - Preserves original CRS in output

        Examples
        --------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>>
        >>> # Hazard data with buffer distances
        >>> data = {
        ...     'ID_hazard': ['h1', 'h2'],
        ...     'buffer_dist_500': [500, 1000],
        ...     'buffer_dist_1000': [1000, 2000],
        ...     'geometry': [Point(-74, 40.7), Point(-74, 40.8)],
        ...     'utm_projection': ['EPSG:32618', 'EPSG:32618']
        ... }
        >>> gdf = gpd.GeoDataFrame(data, crs='EPSG:4326')
        >>>
        >>> # Add buffered geometries
        >>> result = GeometryOperations.add_buffered_geometry_columns(gdf)
        >>> 'buffered_hazard_500' in result.columns
        True
        >>> 'buffered_hazard_1000' in result.columns
        True
        """
        gdf = gdf.copy()

        buffer_cols = [col for col in gdf.columns if col.startswith("buffer_dist")]
        for buffer_col in buffer_cols:
            suffix = buffer_col.replace("buffer_dist", "").strip("_")
            new_col = f"buffered_hazard_{suffix}" if suffix else "buffered_hazard"
            gdf[new_col] = gdf.apply(
                lambda row: GeometryOperations.get_buffered_geometry(row, buffer_col),
                axis=1,
            )

        return gdf

    @staticmethod
    def combine_geometries_by_column(
        gdf: gpd.GeoDataFrame,
        chunk_size: int = 500,
    ) -> gpd.GeoDataFrame:
        """
        Combine all geometries in columns starting with 'buffered_hazard' into a
        single geometry.

        Uses chunking for efficiency with large datasets. Each geometry
        column is processed separately to create a single unioned geometry.

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Input GeoDataFrame
        chunk_size : int, optional
            Number of geometries to process per chunk, by default 500

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame with one row and merged geometry columns

        Examples
        --------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>>
        >>> # Multiple buffered hazards
        >>> data = {
        ...     'ID_hazard': ['h1', 'h2'],
        ...     'buffered_hazard_500': [Point(0, 0).buffer(0.5), Point(2, 2).buffer(0.5)],
        ...     'buffered_hazard_1000': [Point(0, 0).buffer(1.0), Point(2, 2).buffer(1.0)]
        ... }
        >>> gdf = gpd.GeoDataFrame(data, crs='EPSG:4326')
        >>> gdf = gdf.set_geometry('buffered_hazard_500')
        >>>
        >>> # Combine geometries
        >>> combined = GeometryOperations.combine_geometries_by_column(gdf)
        >>> len(combined)
        1
        """
        buffered_cols = [
            col for col in gdf.columns if col.startswith("buffered_hazard")
        ]
        merged_geoms = {}
        for col in buffered_cols:
            geoms = [
                g for g in gdf[col] if g is not None and g.is_valid and not g.is_empty
            ]
            chunks = [
                geoms[i : i + chunk_size] for i in range(0, len(geoms), chunk_size)
            ]
            partial_unions = [unary_union(chunk) for chunk in chunks]
            final_union = unary_union(partial_unions)
            merged_geoms[col] = [final_union]
        merged_geoms["ID_hazard"] = ["merged_geoms"]
        combined_gdf = gpd.GeoDataFrame(
            merged_geoms, geometry=buffered_cols[0], crs=gdf.crs
        )
        return combined_gdf

    @staticmethod
    def get_geometry_intersections(
        hazards_gdf: gpd.GeoDataFrame,
        admin_units_gdf: gpd.GeoDataFrame,
    ) -> gpd.GeoDataFrame:
        """
        Calculate intersections between hazards and administrative units for multiple geometry columns.

        Parameters
        ----------
        hazards_gdf : geopandas.GeoDataFrame
            Hazards GeoDataFrame with ID_hazard column
        admin_units_gdf : geopandas.GeoDataFrame
            Administrative units GeoDataFrame with ID_admin_unit column

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame with intersections for each geometry column

        Examples
        --------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point, Polygon
        >>>
        >>> # Hazard data
        >>> hazards = gpd.GeoDataFrame({
        ...     'ID_hazard': ['h1'],
        ...     'buffered_hazard_500': [Point(0, 0).buffer(0.5)],
        ...     'buffered_hazard_1000': [Point(0, 0).buffer(1.0)]
        ... }, crs='EPSG:4326')
        >>>
        >>> # Admin units
        >>> units = gpd.GeoDataFrame({
        ...     'ID_admin_unit': ['u1'],
        ...     'geometry': [Polygon([(-1, -1), (1, -1), (1, 1), (-1, 1)])]
        ... }, crs='EPSG:4326')
        >>>
        >>> # Get intersections
        >>> intersections = GeometryOperations.get_geometry_intersections(
        ...     hazards, units)
        >>> 'ID_hazard' in intersections.columns
        True
        >>> 'ID_admin_unit' in intersections.columns
        True
        """
        intersections = {}
        for col in [c for c in hazards_gdf.columns if c.startswith("buffered_hazard")]:
            # Select only ID_hazard and the current geometry column
            hazards_subset = hazards_gdf[["ID_hazard", col]].copy()
            hazards_geom = hazards_subset.set_geometry(col, crs=hazards_gdf.crs)
            intersection = gpd.overlay(
                hazards_geom, admin_units_gdf, how="intersection"
            )
            intersection = GeometryValidator.remove_missing_geometries(intersection)
            intersection = GeometryValidator.clean_geometries(intersection)
            intersection = intersection.rename_geometry(col)
            intersection = intersection.set_geometry(col, crs=hazards_gdf.crs)

            intersections[col] = intersection
        intersected_dfs = [
            df for df in intersections.values() if df is not None and not df.empty
        ]

        intersected_hazards = functools.reduce(
            lambda left, right: pd.merge(
                left, right, on=["ID_hazard", "ID_admin_unit"], how="outer"
            ),
            intersected_dfs,
        )
        return intersected_hazards
