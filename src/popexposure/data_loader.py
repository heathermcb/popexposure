"""
Data loading and validation for geospatial population exposure analysis.

This module provides classes for reading and validating geospatial data files
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
        >>> gdf = DataReader.read_geospatial_file("data/hazards.geojson")
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
          - ``ID_hazard``: unique identifier for each hazard (string/object)
          - ``geometry``: spatial geometry (geometry dtype)
          - ``buffer_dist_*``: one or more buffer distance columns (numeric)

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Input GeoDataFrame to validate

        Returns
        -------
        bool
            True if all required columns are present with correct types, False otherwise
        """
        required = ["ID_hazard", "geometry"]
        buffer_cols = [col for col in gdf.columns if col.startswith("buffer_dist")]

        # Check columns exist
        if not (all(col in gdf.columns for col in required) and buffer_cols):
            return False

        # Check types
        try:
            return (
                gdf["ID_hazard"].dtype == "object"
                and gdf.geometry.dtype.name == "geometry"
                and all(gdf[col].dtype.kind in "biufc" for col in buffer_cols)
            )
        except (KeyError, AttributeError):
            return False

    @staticmethod
    def validate_admin_unit_columns(gdf: gpd.GeoDataFrame) -> bool:
        """
        Validate that GeoDataFrame has required columns for admin unit analysis.

        Expected columns:
          - ``ID_admin_unit``: unique identifier for each admin unit (string/object)
          - ``geometry``: admin unit geometry (geometry dtype)

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Input GeoDataFrame to validate

        Returns
        -------
        bool
            True if all required columns are present with correct types, False otherwise
        """
        required = ["ID_admin_unit", "geometry"]

        # Check columns exist
        if not all(col in gdf.columns for col in required):
            return False

        # Check types
        try:
            return (
                gdf["ID_admin_unit"].dtype == "object"
                and gdf.geometry.dtype.name == "geometry"
            )
        except (KeyError, AttributeError):
            return False
