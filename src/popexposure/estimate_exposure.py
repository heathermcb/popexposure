"""
Population exposure estimation for environmental hazards.

Main interface for calculating populations exposed to environmental hazards
using geospatial analysis and gridded population data.
"""

import geopandas as gpd
import pandas as pd

from .data_loader import DataReader as dr
from .geometry_validator import GeometryValidator as gv
from .geometry_operations import GeometryOperations as go
from .raster_extraction import RasterExtractor as re


class PopEstimator:
    """
    Estimate population exposure to environmental hazards using geospatial analysis.

    PopEstimator provides a complete workflow for calculating how many people live
    within specified buffer distances of environmental hazards (e.g., wildfires,
    oil wells, toxic sites) using gridded population data. The class handles data
    loading, geometry processing, buffering operations, and raster value extraction
    to produce exposure estimates.

    Key Features
    ------------
    - **Flexible hazard data**: Works with point, line, polygon, multipolygon, or geometry collection hazards
    - **Multiple buffer distances**: Calculate exposure at different distances simultaneously
    - **Administrative breakdowns**: Get exposure counts by census tracts, ZIP codes, etc.
    - **Hazard-specific or combined estimates**: Choose individual hazard impacts or cumulative exposure (see est_exposed_pop)
    - **Automatic geometry processing**: Handles CRS transformations, invalid geometries, and projections seamlessly
    - **Partial pixel extraction**: Uses area-weighted raster sampling for accurate population counts

    Workflow
    --------
    1. **Load and clean data** with :meth:`prep_data`
    2. **Calculate exposure** with :meth:`est_exposed_pop`
    3. **Get total administrative unit populations** with :meth:`est_pop` (optional)

    Parameters
    ----------
    No parameters required for initialization.

    Attributes
    ----------
    hazard_data : geopandas.GeoDataFrame or None
        Processed hazard data with buffered geometries (set by prep_data)
    admin_units : geopandas.GeoDataFrame or None
        Administrative unit geometries (set by prep_data)
    population : pandas.DataFrame or None
        Total population by administrative unit (set by est_pop)

    Examples
    --------
    Basic exposure analysis:

    >>> import popexposure
    >>>
    >>> # Initialize estimator
    >>> estimator = popexposure.PopEstimator()
    >>>
    >>> # Load hazard data (e.g., oil wells with 500m and 1000m buffers)
    >>> hazards = estimator.prep_data("oil_wells.geojson", "hazard")
    >>>
    >>> # Calculate population exposure
    >>> exposure = estimator.est_exposed_pop(
    ...     pop_path="population_raster.tif",
    ...     hazard_specific=True  # Individual hazard estimates
    ... )
    >>> print(exposure.head())
       ID_hazard  exposed_500  exposed_1000
    0      well_1         1234          3456
    1      well_2          856          2134

    Administrative breakdown:

    >>> # Load administrative units (e.g., census tracts)
    >>> admin_units = estimator.prep_data("census_tracts.geojson", "admin_unit")
    >>>
    >>> # Get exposure by census tract
    >>> exposure_by_tract = estimator.est_exposed_pop(
    ...     pop_path="population_raster.tif",
    ...     hazard_specific=False,  # Combined hazard exposure
    ...     admin_units=admin_units
    ... )
    >>>
    >>> # Get total population for percentage calculations
    >>> total_pop = estimator.est_pop("population_raster.tif", admin_units)
    >>>
    >>> # Calculate exposure percentages
    >>> exposure_pct = exposure_by_tract.merge(total_pop, on="ID_admin_unit")
    >>> exposure_pct["pct_exposed_500"] = (
    ...     exposure_pct["exposed_500"] / exposure_pct["population"] * 100
    ... )

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
    prep_data : Load and preprocess geospatial data
    est_exposed_pop : Calculate population exposure to hazards
    est_pop : Calculate total population in administrative units
    """

    def __init__(self):
        """
        Initialize the PopEstimator class, used to find populations exposed to
        environmental hazards.
        Init with empty attributes for hazard and admin unit data.
        """
        self.hazard_data = None
        self.admin_units = None
        self.population = None
        self.exposed = None

    def prep_data(self, path_to_data: str, geo_type: str) -> gpd.GeoDataFrame:
        """
        Reads, cleans, and preprocesses geospatial data for exposure analysis.

        This method loads a geospatial file (GeoJSON or GeoParquet) containing
        either hazard data (e.g., wildfire burn zones, oil wells)
        or additional administrative geographies (e.g., ZCTAs, census tracts,
        referred to here as admin units). It makes all geometries valid,
        removes empty geometries, and, for hazard data, generates buffered
        geometries for one or more user-specified buffer distances.
        Buffering is performed in the best Universal Transverse Mercator (UTM)
        projection based on each geometry's centroid latitude and longitude.

        Parameters
        ----------
        path_to_data : str
            Path to a geospatial data file (.geojson or .parquet). The file must
            contain either hazard data or administrative geography data, as
            specified by ``geo_type``.
            Data must have any coordinate reference system.
            - Hazard data must contain a string column ``"ID_hazard"`` with
            unique hazard IDs, a geometry column ``"geometry"``, and one or more
            numeric columns starting with ``"buffer_dist"`` with unique suffixes
            (e.g., ``"buffer_dist_main"``, ``"buffer_dist_1000"``) specifying
            buffer distances in meters. Buffer distances may be 0 or different
            for each hazard.
            - For administrative unit data, the file must contain a string column
            ``"ID_admin_unit"`` with unique admin unit IDs and a geometry column
            ``"geometry"``.
        geo_type : str
            A string indicating the type of data to process. Must be either
            ``"hazard"`` for environmental hazard data or ``"admin_unit"`` for
            administrative geography data.


        Returns
        -------
        geopandas.GeoDataFrame or None
            A GeoDataFrame with cleaned and valid geometries.
            - If hazard data was passed, the output contains a column
            ``"ID_hazard"`` matching the input data, and one or more
            ``"buffered_hazard"`` geometry columns, with suffixes matching the
            passed ``buffer_dist`` columns (e.g., ``"buffered_hazard_main"``,
            ``"buffered_hazard_1000"``).
            - If admin unit data was passed, the output contains columns
            ``"ID_admin_unit"`` matching the input data and ``"geometry"``.
            - Empty geometries are removed.
            - If the input file is empty or contains no valid geometries, the
            function returns None.
        """
        shp_df = dr.read_geospatial_file(path_to_data)
        if shp_df.empty:
            return None

        shp_df = gv.remove_missing_geometries(shp_df)
        shp_df = gv.clean_geometries(shp_df)
        shp_df = gv.reproject_to_wgs84(shp_df)

        if geo_type == "hazard":
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
            self.hazard_data = buffered_hazards
            return buffered_hazards

        elif geo_type == "admin_unit":
            self.admin_units = shp_df
            return shp_df

        else:
            raise ValueError("geo_type must be 'hazard' or 'admin_unit'")

    def est_exposed_pop(
        self,
        pop_path: str,
        hazard_specific: bool,
        hazards: gpd.GeoDataFrame = None,
        admin_units: gpd.GeoDataFrame = None,
        stat: str = "sum",
    ) -> pd.DataFrame:
        """
        Estimate the number of people living within a buffer distance of
        environmental hazard(s) using a gridded population raster.

        This function calculates the sum of raster values within buffered hazard
        geometries, or within the intersection of buffered hazard geometries and
        additional administrative geographies, to find the population exposed to
        hazards. Users can choose to estimate either (a) hazard-specific counts
        (the number of people exposed to each unique buffered hazard in the set),
        or (b) a cumulative count (the number of unique people exposed to any
        of the input buffered hazards, avoiding double counting). Either
        estimate can be broken down by additional administrative geographies
        such as ZCTAs. Users must supply at least one buffered hazard column,
        but may supply additional buffered hazard columns to create estimates
        of exposure for different buffer distances.

        Parameters
        ----------
        pop_path : str
            Path to a gridded population raster file, in TIFF format. The raster
            must have any coordinate reference system.
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
        admin_units : geopandas.GeoDataFrame, optional
            An optional GeoDataFrame of additional administrative geographies,
            containing a string column called ``ID_admin_unit`` and a geometry
            column called ``geometry``.
        stat : str, default "sum"
            Statistic to calculate from raster values. Options:
            - "sum": Total population within geometry (default)
            - "mean": Average raster value within geometry

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

            The number of rows in the output DataFrame depends on the function
            arguments:
            - If ``hazard_specific`` is True, the DataFrame contains one row per
            hazard or per hazard-admin unit pair, if admin units are provided.
            - If ``hazard_specific`` is False, the DataFrame contains a single
            row or one row per admin unit, if admin units are provided, with
            each ``exposed`` column representing the total population in the
            union of all buffered hazard geometries in that buffered hazard column.

        Notes
        -----
        There are four ways to use this function:

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
           merged into a single geometry, and the function calculates the total
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

        if hazards is None:
            hazards = self.hazard_data
        if admin_units is None:
            admin_units = self.admin_units
        if hazards is None:
            return None

        if admin_units is None:
            if not hazard_specific:
                hazards = go.combine_geometries_by_column(hazards)
            exposed = re.mask_raster_partial_pixel(hazards, pop_path, stat=stat)
            self.exposed = exposed
            return exposed

        else:
            if not hazard_specific:
                hazards = go.combine_geometries_by_column(hazards)
            intersected_hazards = go.get_geometry_intersections(
                hazards_gdf=hazards, admin_units_gdf=admin_units
            )
            exposed = re.mask_raster_partial_pixel(
                intersected_hazards, raster_path=pop_path, stat=stat
            )
            self.exposed = exposed
            return exposed

    def est_pop(
        self, pop_path: str, admin_units: str, stat: str = "sum"
    ) -> pd.DataFrame:
        """
        Estimate the total population residing within administrative geographies
        using a gridded population raster.

        This function estimates the total population residing within administrative
        geographies (e.g., ZCTAs, census tracts) according to a provided gridded
        population raster. This method is meant to be used with the same population
        raster as ``est_exposed_pop`` to provide denominators for the total population
        in each administrative geography, allowing the user to compute the
        percentage of people exposed to hazards in each admin unit. ``est_pop``
        calculates the sum of raster values within the boundaries of each
        administrative geography geometry provided.

        Parameters
        ----------
        pop_path : str
            Path to a gridded population raster file, in TIFF format. The raster
            must have any coordinate reference system.
        admin_units : geopandas.GeoDataFrame
            GeoDataFrame containing administrative geography geometries. Must
            include a string column called ``ID_admin_unit`` with unique admin
            unit IDs and a geometry column called ``geometry``.
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
        residing = re.mask_raster_partial_pixel(
            admin_units, raster_path=pop_path, stat=stat
        )
        residing = residing.rename(
            columns=lambda c: c.replace("exposedgeometry", "population")
        )
        self.population = residing
        return residing
