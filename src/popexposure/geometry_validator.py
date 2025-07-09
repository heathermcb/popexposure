"""
Geometry validation and coordinate reference system utilities.

This module provides classes for basic geometry validation, cleaning,
and CRS operations used in population exposure calculations.
"""

import geopandas as gpd
from shapely.validation import make_valid
from shapely.geometry import Point


class GeometryValidator:
    """Handles geometry validation, cleaning, and CRS operations."""

    @staticmethod
    def remove_missing_geometries(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Remove rows with null or empty geometries.

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Input GeoDataFrame

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame with missing geometries removed

        Examples
        --------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>> import numpy as np
        >>>
        >>> # Data with missing geometries
        >>> data = {
        ...     'ID': ['a', 'b', 'c'],
        ...     'geometry': [Point(0, 0), None, Point(1, 1)]
        ... }
        >>> gdf = gpd.GeoDataFrame(data)
        >>> cleaned = GeometryValidator.remove_missing_geometries(gdf)
        >>> len(cleaned)
        2
        """
        valid_mask = gdf["geometry"].notnull() & ~gdf["geometry"].is_empty

        if not valid_mask.all():
            id_col = next(col for col in gdf.columns if "ID" in col)
            dropped_ids = gdf.loc[~valid_mask, id_col].tolist()
            print(
                f"Warning: {len(dropped_ids)} geometries with null/empty values were dropped. IDs: {dropped_ids}"
            )

        return gdf[valid_mask].copy().reset_index(drop=True)

    @staticmethod
    def clean_geometries(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Clean geometries by making them valid using buffer(0) and make_valid().

        This method applies both buffer(0) and make_valid() to fix common
        geometry issues including self-intersections, duplicate vertices,
        and topology problems. Geometries that cannot be made valid are
        removed and a warning is printed.

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Input GeoDataFrame

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame with cleaned geometries

        Examples
        --------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Polygon
        >>>
        >>> # Self-intersecting polygon (bowtie)
        >>> coords = [(0, 0), (2, 2), (2, 0), (0, 2), (0, 0)]
        >>> invalid_poly = Polygon(coords)
        >>> gdf = gpd.GeoDataFrame({'ID': ['valid_geom'], 'geometry': [invalid_poly]})
        >>>
        >>> # Clean geometries
        >>> cleaned_gdf = GeometryValidator.clean_geometries(gdf)
        >>> cleaned_gdf.geometry.is_valid.all()
        True
        """
        gdf = gdf.copy()

        # Apply make_valid to handle any remaining issues
        gdf["geometry"] = gdf["geometry"].apply(make_valid)

        # Check for geometries that are still invalid or became empty
        valid_mask = (
            gdf["geometry"].notnull()
            & gdf["geometry"].is_valid
            & ~gdf["geometry"].is_empty
        )

        # Print warning for dropped geometries
        if not valid_mask.all():
            id_col = next(col for col in gdf.columns if "ID" in col)
            dropped_ids = gdf.loc[~valid_mask, id_col].tolist()
            print(
                f"Warning: {len(dropped_ids)} geometries could not be made valid and were dropped. IDs: {dropped_ids}"
            )

        return gdf[valid_mask].copy().reset_index(drop=True)

    @staticmethod
    def reproject_to_wgs84(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Reproject all geometries to WGS84 (EPSG:4326).

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Input GeoDataFrame

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame with geometries reprojected to WGS84

        Examples
        --------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>>
        >>> # Data in Web Mercator projection (common for web maps)
        >>> gdf = gpd.GeoDataFrame(
        ...     {'geometry': [Point(-8238310, 4969803)]},
        ...     crs='EPSG:3857'  # Web Mercator
        ... )
        >>>
        >>> # Reproject to WGS84
        >>> wgs84_gdf = GeometryValidator.reproject_to_wgs84(gdf)
        >>> str(wgs84_gdf.crs)
        'EPSG:4326'
        """
        if gdf.crs != "EPSG:4326":
            return gdf.to_crs("EPSG:4326")
        return gdf

    @staticmethod
    def get_best_utm_projection(lat: float, lon: float) -> str:
        """
        Get the best UTM projection EPSG code for given coordinates.

        Parameters
        ----------
        lat : float
            Latitude in decimal degrees
        lon : float
            Longitude in decimal degrees

        Returns
        -------
        str
            EPSG code string for the best UTM projection

        Examples
        --------
        >>> # New York City coordinates
        >>> utm_code = GeometryValidator.get_best_utm_projection(40.7128, -74.0060)
        >>> utm_code
        'EPSG:32618'

        >>> # Sydney coordinates
        >>> utm_code = GeometryValidator.get_best_utm_projection(-33.8688, 151.2093)
        >>> utm_code
        'EPSG:32756'
        """
        zone_number = (lon + 180) // 6 + 1
        hemisphere = 326 if lat >= 0 else 327
        epsg_code = hemisphere * 100 + zone_number
        return f"EPSG:{int(epsg_code)}"

    @staticmethod
    def add_utm_projection_column(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Add UTM projection column based on geometry centroids.

        Calculates centroid coordinates and determines the optimal UTM
        projection for each geometry.

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Input GeoDataFrame with geometries in WGS84

        Returns
        -------
        geopandas.GeoDataFrame
            GeoDataFrame with added column: utm_projection

        Examples
        --------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>>
        >>> # Geometries in different locations
        >>> data = {
        ...     'ID_hazard': ['h1', 'h2'],
        ...     'geometry': [Point(-74, 40.7), Point(151, -33.9)]  # NYC, Sydney
        ... }
        >>> gdf = gpd.GeoDataFrame(data, crs='EPSG:4326')
        >>>
        >>> # Add UTM projections
        >>> result = GeometryValidator.add_utm_projection_column(gdf)
        >>> 'utm_projection' in result.columns
        True
        >>> result['utm_projection'].iloc[0]
        'EPSG:32618'
        >>> result['utm_projection'].iloc[1]
        'EPSG:32756'
        """
        gdf = gdf.copy()

        # reproject to molleweide for centroid calc
        gdf = gdf.to_crs("ESRI:54009")

        # Get centroid coordinates in Mollweide (meters)
        mollweide_centroids_x = gdf.centroid.x
        mollweide_centroids_y = gdf.centroid.y

        # Create temporary GeoDataFrame with centroid points in Mollweide
        centroid_points = gpd.GeoDataFrame(
            geometry=[
                Point(x, y)
                for x, y in zip(mollweide_centroids_x, mollweide_centroids_y)
            ],
            crs="ESRI:54009",  # Still in Mollweide
        )

        # Reproject the centroid points to WGS84
        centroid_points_wgs84 = centroid_points.to_crs("EPSG:4326")

        # Extract the WGS84 coordinates for UTM zone calculation
        gdf["centroid_lon"] = centroid_points_wgs84.geometry.x
        gdf["centroid_lat"] = centroid_points_wgs84.geometry.y

        # reproject main geometries back to wgs84
        gdf = gdf.to_crs("EPSG:4326")

        # Get UTM projection for each geometry using proper lat/lon coordinates
        gdf["utm_projection"] = gdf.apply(
            lambda row: GeometryValidator.get_best_utm_projection(
                lat=row["centroid_lat"], lon=row["centroid_lon"]
            ),
            axis=1,
        )

        # Remove temporary centroid coordinate columns
        gdf = gdf.drop(columns=["centroid_lon", "centroid_lat"])

        return gdf
