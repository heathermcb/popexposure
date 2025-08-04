"""
Population exposure estimation for environmental hazards.

Main interface for calculating populations exposed to environmental hazards
using geospatial analysis and gridded population data.
"""

from pathlib import Path
import geopandas as gpd
import pandas as pd
from typing import Literal

from .utils import reader as rdr
from .utils import geom_validator as gv
from .utils import geom_ops as go
from .utils.mask_raster_partial_pixel import mask_raster_partial_pixel


class PopEstimator:
    """
    Estimate population exposure to environmental hazards using geospatial analysis.

    PopEstimator provides a complete workflow for calculating how many people live
    within specified buffer distances of environmental hazards (e.g., wildfires,
    oil wells, toxic sites) using gridded population data. The class handles data
    loading, geometry processing, buffering operations, and raster value extraction
    to produce exposure estimates.

    Parameters
    ----------
    One parameter is required for initialization, another is optional.


    Attributes
    ----------
    hazard_data : geopandas.GeoDataFrame or None
        Processed hazard data with buffered geometries (set by prep_data)
    admin_units : geopandas.GeoDataFrame or None
        Administrative unit geometries (set by prep_data)
    population : pandas.DataFrame or None
        Total population by administrative unit (set by est_total_pop)

    Key Features
    ------------
    - **Flexible hazard data**: Works with point, line, polygon, multipolygon, or geometry collection hazards.
    - **Multiple buffer distances**: Calculate exposure at different distances simultaneously.
    - **Administrative breakdowns**: Get exposure counts by census tracts, ZIP codes, etc.
    - **Hazard-specific or combined estimates**: Choose individual hazard impacts or cumulative exposure (see est_exposed_pop).
    - **Automatic geometry processing**: Handles CRS transformations, invalid geometries, and projections seamlessly.
    - **Partial pixel extraction**: Uses area-weighted raster sampling for accurate population counts.

    Workflow
    --------
    1. **Construct an estimator containing population and administrative data.
    2. **Calculate exposure** with :meth:`est_exposed_pop`.
    3. **Get total administrative unit populations** with :meth:`est_total_pop` (optional).

    Examples
    --------
    Basic exposure analysis without admin data:

    >>> import popexposure
    >>>
    >>> # Initialize with only population raster (no admin data)
    >>> estimator = PopEstimator(
    ...     pop_data="data/population.tif"
    ... )
    >>> # Estimate population exposure to hazards (e.g., wildfires)
    >>> exposure = estimator.est_exposed_pop(
    ...     hazard_data="data/wildfire_perimeters.geojson",
    ...     hazard_specific=True
    ... )
    >>> print(exposure.head())
    ID_hazard  exposed_500  exposed_1000
    0   fire_001        1234          2345
    1   fire_002         567           890


    Exposure analysis with admin data:

    >>> import popexposure
    >>> # Initialize with population raster and admin boundaries
    >>> estimator = PopEstimator(
    ...     pop_data="data/population.tif",
    ...     admin_data="data/admin_units.geojson"
    ... )
    >>> # Estimate population exposure to hazards (e.g., wildfires)
    >>> exposure = estimator.est_exposed_pop(
    ...     hazard_data="data/wildfire_perimeters.geojson",
    ...     hazard_specific=True
    ... )
    >>> print(exposure.head())
    ID_hazard  ID_admin_unit  exposed_500  exposed_1000
    0   fire_001         06001        1234          2345
    1   fire_002         06013         567           890

    >>> # Estimate total population in each admin unit
    >>> total_pop = estimator.est_total_pop()
    >>> print(total_pop.head())
    ID_admin_unit  population
    0         06001      100000
    1         06013      150000

    Notes
    -----
    **Data Requirements:**

    - **Hazard data**: Must contain ``ID_hazard`` column ``buffer_dist_*`` columns, and ``geometry`` column
    - **Admin units**: Must contain ``ID_admin_unit`` column and ``geometry`` column
    - **Population raster**: Any CRS supported

    **Buffer Distance Naming:**

    - Column ``buffer_dist_500`` creates ``buffered_hazard_500`` and ``exposed_500``
    - Column ``buffer_dist_main`` creates ``buffered_hazard_main`` and ``exposed_main``
    - Distances are in meters and can vary by hazard

    **Coordinate Reference Systems:**

    - Input data can use any CRS
    - Buffering uses optimal UTM projections for accuracy
    - Population raster CRS is automatically handled

    See Also
    --------
    est_exposed_pop : Calculate population exposure to hazards
    est_total_pop : Calculate total population in administrative units
    """

    def __init__(
        self,
        pop_data: str | Path,
        admin_data: str | Path | gpd.GeoDataFrame | None = None,
    ):
        """
        Initialize the PopEstimator class, used to find populations exposed to
        environmental hazards, with a population raster and optional admin units.
        """
        self.pop_data = pop_data
        if admin_data is not None:
            if isinstance(admin_data, gpd.GeoDataFrame) and self._is_admin_data_prepped(
                admin_data
            ):
                self.admin_data = admin_data.copy()
            else:
                self.admin_data = self._process_admin_data(admin_data)
        else:
            self.admin_data = None

    def _is_admin_data_prepped(self, admin_gdf: gpd.GeoDataFrame) -> bool:
        """
        Check if admin data is already prepped and ready to use.

        Prepped admin data should have:
        - An ID column (column name containing 'ID')
        - No null/empty geometries
        - Clean geometries (valid according to Shapely)
        - WGS84 projection (EPSG:4326)

        Parameters
        ----------
        admin_gdf : geopandas.GeoDataFrame
            Admin data to check

        Returns
        -------
        bool
            True if data is prepped, False otherwise
        """
        # Check for ID column
        id_cols = [col for col in admin_gdf.columns if "ID" in col]
        if not id_cols:
            return False

        # Check for valid geometries (no null/empty)
        if admin_gdf.geometry.isnull().any() or admin_gdf.geometry.is_empty.any():
            return False

        # Check for WGS84 projection
        if admin_gdf.crs is None or admin_gdf.crs.to_epsg() != 4326:
            return False

        # Check for valid geometries (Shapely validity)
        if not admin_gdf.geometry.is_valid.all():
            return False

        return True

    def _is_hazard_data_prepped(self, hazard_gdf: gpd.GeoDataFrame) -> bool:
        """
        Check if hazard data is already prepped and ready to use.

        Prepped hazard data should have:
        - An ID_hazard column
        - No null/empty geometries
        - Clean geometries (valid according to Shapely)
        - WGS84 projection (EPSG:4326)
        - One or more buffered_hazard_* columns
        - Geometry set to a buffered_hazard column

        Parameters
        ----------
        hazard_gdf : geopandas.GeoDataFrame
            Hazard data to check

        Returns
        -------
        bool
            True if data is prepped, False otherwise
        """
        # Check for ID_hazard column
        if "ID_hazard" not in hazard_gdf.columns:
            return False

        # Check for buffered_hazard columns
        buffered_cols = [
            col for col in hazard_gdf.columns if col.startswith("buffered_hazard")
        ]
        if not buffered_cols:
            return False

        # Check for valid geometries (no null/empty)
        if hazard_gdf.geometry.isnull().any() or hazard_gdf.geometry.is_empty.any():
            return False

        # Check for WGS84 projection
        if hazard_gdf.crs is None or hazard_gdf.crs.to_epsg() != 4326:
            return False

        # Check for valid geometries (Shapely validity)
        if not hazard_gdf.geometry.is_valid.all():
            return False

        # Check that geometry is set to a buffered_hazard column
        geom_col = hazard_gdf.geometry.name
        if not geom_col.startswith("buffered_hazard"):
            return False

        return True

    def _process_admin_data(self, data: str | Path | gpd.GeoDataFrame):
        shp_df = (
            rdr.read_geospatial_file(data)
            if isinstance(data, (str, Path))
            else data.copy()
        )
        shp_df = gv.remove_missing_geometries(shp_df)
        shp_df = gv.clean_geometries(shp_df)
        shp_df = gv.reproject_to_wgs84(shp_df)
        return shp_df

    def _process_hazard_data(self, data: str | Path | gpd.GeoDataFrame):
        shp_df = (
            rdr.read_geospatial_file(data)
            if isinstance(data, (str, Path))
            else data.copy()
        )

        shp_df = gv.remove_missing_geometries(shp_df)
        shp_df = gv.clean_geometries(shp_df)
        shp_df = gv.reproject_to_wgs84(shp_df)

        shp_df = gv.add_utm_projection_column(shp_df)
        shp_df = go.add_buffered_geometry_columns(shp_df)
        buffered_cols = [
            col for col in shp_df.columns if col.startswith("buffered_hazard")
        ]
        cols = ["ID_hazard"] + buffered_cols
        buffered_hazards = shp_df[cols]
        buffered_hazards = buffered_hazards.set_geometry(
            buffered_cols[0], crs="EPSG:4326"
        )
        return buffered_hazards

    def est_exposed_pop(
        self,
        hazard_data: str | Path | gpd.GeoDataFrame,
        hazard_specific: bool = True,
        stat: Literal["sum", "mean"] = "sum",
    ) -> pd.DataFrame:
        """
        Estimate the number of people living within a buffer distance of
        environmental hazard(s) using a gridded population raster.

        This method calculates the population exposed to hazards by summing
        raster values within buffered hazard geometries, or within the
        intersection of these buffers and administrative geographies (if
        provided to the class). Users can choose between hazard-specific
        counts (population exposed to each individual hazard) or cumulative
        counts (population exposed to any hazard, without double counting).
        Exposure can be estimated for multiple buffer distances simultaneously,
        as specified by the buffered hazard columns in the input data. If
        administrative units were supplied when initializing the class, results
        are further broken down by these geographies (e.g., census tracts or
        ZIP codes). At least one buffered hazard column must be present in the
        hazard data; additional columns allow for exposure estimates at multiple
        distances.

        Parameters
        ----------

        hazard_specific : bool
            If True, exposure is calculated for each hazard individually
            (hazard-specific estimates). If False, geometries are combined before
            exposure is calculated, producing a single cumulative estimate.
        hazards : geopandas.GeoDataFrame
            A GeoDataFrame with a coordinate reference system containing a
            string column called ``ID_hazard`` with unique hazard IDs, and one
            or more geometry columns starting with ``buffered_hazard``
            containing buffered hazard geometries. ``buffered_hazard`` columns
            must each have a unique suffix (e.g., ``buffered_hazard_10``,
            ``buffered_hazard_100``, ``buffered_hazard_1000``).
        stat : str, default "sum"
            Statistic to calculate from raster values. Options:
            - "sum": Total population within geometry (default)
            - "mean": Average raster value/population value within geometry

        Returns
        -------
        pandas.DataFrame
            A DataFrame with the following columns:
            - ``ID_hazard``: Always included.
            - ``ID_admin_unit``: Included only if admin units were provided.
            - One or more ``exposed`` columns: Each corresponds to a buffered
            hazard column (e.g., if the input had columns ``buffered_hazard_10``,
            ``buffered_hazard_100``, and ``buffered_hazard_1000``, the output
            will have ``exposed_10``, ``exposed_100``, and ``exposed_1000``).
            Each ``exposed`` column contains the statistic (sum or mean) of raster
            values (population) within the relevant buffered hazard geometry or
            buffered hazard geometry and admin unit intersection.

            The number of rows in the output DataFrame depends on the method
            arguments:
            - If ``hazard_specific`` is True, the DataFrame contains one row per
            hazard or per hazard-admin unit pair, if admin units are provided.
            - If ``hazard_specific`` is False, the DataFrame contains a single
            row or one row per admin unit, if admin units are provided, with
            each ``exposed`` column representing the total population in the
            union of all buffered hazard geometries in that buffered hazard column.

        Notes
        -----
        There are four ways to use this method:

        1. Hazard-specific exposure, no additional administrative geographies
        (``hazard_specific=True, admin_units=None``):
           Calculates the exposed population for each buffered hazard geometry.
           Returns a DataFrame with one row per hazard and one ``exposed``
           column per buffered hazard column. If people lived within the buffer
           distance of more than one hazard, they are included in the exposure
           counts for each hazard they are near.

        2. Combined hazards, no additional administrative geographies
        (``hazard_specific=False, admin_units=None``):
           All buffered hazard geometries in each buffered hazard column are
           merged into a single geometry, and the method calculates the total
           exposed population for the union of those buffered hazards. Returns a
           DataFrame with a single row and one ``exposed`` column for each
           buffered hazard column. If people were close to more than one hazard
           in the hazard set, they are counted once.

        3. Hazard-specific exposure within admin units
        (``hazard_specific=True, admin_units`` provided):
           Calculates the exposed population for each intersection of each
           buffered hazard geometry and each admin unit. Returns a DataFrame
           with one row per buffered hazard-admin unit pair and one ``exposed``
           column per buffered hazard column. If people lived within the buffer
           distance of more than one hazard, they are included in the exposure
           counts for their admin unit-hazard combination for each hazard
           they are near.

        4. Combined hazards within admin units
        (``hazard_specific=False, admin_units`` provided):
           All buffered hazard geometries in the same column are merged into a
           single geometry. Calculates the exposed population for the
           intersection of each buffered hazard combined geometry with each admin
           unit. Returns a DataFrame with one row per admin unit and one
           ``exposed`` column per buffered hazard column. If people were close
           to more than one hazard in the hazard set, they are counted once.
        """
        if isinstance(hazard_data, gpd.GeoDataFrame) and self._is_hazard_data_prepped(
            hazard_data
        ):
            hazard_data = hazard_data.copy()
        else:
            hazard_data = self._process_hazard_data(hazard_data)

        if not hazard_specific:
            hazard_data = go.combine_geometries_by_column(hazard_data)

        if self.admin_data is not None:
            hazard_data = go.get_geometry_intersections(
                hazards_gdf=hazard_data, admin_units_gdf=self.admin_data
            )

        exposed = mask_raster_partial_pixel(
            hazard_data, raster_path=self.pop_data, stat=stat
        )

        return exposed

    def est_total_pop(self, stat: Literal["sum", "mean"] = "sum") -> pd.DataFrame:
        """
        Estimate the total population residing within administrative geographies
        using a gridded population raster.

        This method estimates the total population residing within administrative
        geographies (e.g., ZCTAs, census tracts) according to a provided gridded
        population raster. This method is meant to be used with the same population
        raster as ``est_exposed_pop`` to provide denominators for the total population
        in each administrative geography, allowing the user to compute the
        percentage of people exposed to hazards in each admin unit. ``est_total_pop``
        calculates the sum of raster values within the boundaries of each
        administrative geography geometry provided.

        Parameters
        ----------
        stat : str, default "sum"
            Statistic to calculate from raster values. Options:
            - "sum": Total population within geometry (default)
            - "mean": Average raster value within geometry

        Returns
        -------
        pandas.DataFrame
            DataFrame with an ``ID_admin_unit`` column matching the input and a
            ``population`` column, where each value is the specified statistic
            (sum or mean) of raster values within the corresponding admin unit geometry.
        """
        residing = mask_raster_partial_pixel(
            self.admin_data, raster_path=self.pop_data, stat=stat
        )
        residing = residing.rename(
            columns=lambda c: c.replace("exposedgeometry", "population")
        )
        return residing


__all__ = ["PopEstimator"]
