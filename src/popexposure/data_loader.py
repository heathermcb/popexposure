"""
Data reading utilities for population exposure analysis.

This module provides classes for reading geospatial data files
used in population exposure calculations.
"""

from pathlib import Path
import geopandas as gpd


class DataReader:
    """Handles reading and validating geospatial data files."""

    @staticmethod
    def read_geospatial_file(path: str) -> gpd.GeoDataFrame:
        """
        Read geospatial data from file.

        Parameters
        ----------
        path : str
            Path to geospatial file (.geojson or .parquet)

        Returns
        -------
        geopandas.GeoDataFrame
            Loaded geospatial data

        Raises
        ------
        FileNotFoundError
            If file type is not supported

        Examples
        --------
        >>> reader = DataReader()
        >>> gdf = reader.read_geospatial_file("data/hazards.geojson")
        >>> print(gdf.shape)
        (100, 3)
        """
        path = Path(path)
        if path.suffix == ".geojson":
            return gpd.read_file(path)
        elif path.suffix == ".parquet":
            return gpd.read_parquet(path)
        else:
            raise FileNotFoundError(f"Unsupported file type: {path}")

    @staticmethod
    def validate_hazard_columns(gdf: gpd.GeoDataFrame) -> bool:
        """
        Validate that GeoDataFrame has required columns for hazard analysis.

        Expected columns:
          - ``ID_hazard``: unique identifier for each hazard
          - ``geometry``: spatial geometry
          - ``buffer_dist_*``: one or more buffer distance columns

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Input GeoDataFrame to validate

        Returns
        -------
        bool
            True if all required columns are present, False otherwise

        Examples
        --------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>>
        >>> # Valid hazard data
        >>> data = {
        ...     'ID_hazard': ['h1', 'h2'],
        ...     'buffer_dist_500': [500, 1000],
        ...     'geometry': [Point(0, 0), Point(1, 1)]
        ... }
        >>> gdf = gpd.GeoDataFrame(data)
        >>> DataReader.validate_hazard_columns(gdf)
        True

        >>> # Invalid hazard data (missing buffer_dist column)
        >>> data = {
        ...     'ID_hazard': ['h1', 'h2'],
        ...     'geometry': [Point(0, 0), Point(1, 1)]
        ... }
        >>> gdf = gpd.GeoDataFrame(data)
        >>> DataReader.validate_hazard_columns(gdf)
        False
        """
        required = ["ID_hazard", "geometry"]
        buffer_cols = [col for col in gdf.columns if col.startswith("buffer_dist")]
        return all(col in gdf.columns for col in required) and len(buffer_cols) > 0

    @staticmethod
    def validate_admin_unit_columns(gdf: gpd.GeoDataFrame) -> bool:
        """
        Validate that GeoDataFrame has required columns for admin unit analysis.

        Expected columns:
          - ``ID_admin_unit``: unique identifier for each admin unit
          - ``geometry``: admin unit geometry

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Input GeoDataFrame to validate

        Returns
        -------
        bool
            True if all required columns are present, False otherwise

        Examples
        --------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Polygon
        >>>
        >>> # Valid spatial unit data
        >>> data = {
        ...     'ID_admin_unit': ['unit1', 'unit2'],
        ...     'geometry': [
        ...         Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
        ...         Polygon([(1, 0), (2, 0), (2, 1), (1, 1)])
        ...     ]
        ... }
        >>> gdf = gpd.GeoDataFrame(data)
        >>> DataReader.validate_admin_unit_columns(gdf)
        True

        >>> # Invalid admin unit data (missing ID_admin_unit)
        >>> data = {
        ...     'name': ['unit1', 'unit2'],
        ...     'geometry': [
        ...         Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
        ...         Polygon([(1, 0), (2, 0), (2, 1), (1, 1)])
        ...     ]
        ... }
        >>> gdf = gpd.GeoDataFrame(data)
        >>> DataReader.validate_admin_unit_columns(gdf)
        False
        """
        required = ["ID_admin_unit", "geometry"]
        return all(col in gdf.columns for col in required)
