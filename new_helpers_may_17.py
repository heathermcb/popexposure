import geopandas as gpd
import pandas as pd
import numpy as np
import pyarrow.parquet as pq
import rasterio
import rasterio.windows
from tqdm import tqdm
from shapely.validation import make_valid
from shapely.ops import unary_union, transform
from shapely import wkt
from pathlib import Path
from exactextract import exact_extract
import warnings
import functools
from osgeo import gdal
import pyproj

# Suppress warnings about centroid crs but raise exceptions
gdal.UseExceptions()
warnings.filterwarnings("ignore")


class PopEstimator:

    def __init__(self):
        """
        Initialize the PopEstimator class, used to find populations exposed to
        environmental hazards.
        Init with empty attributes for hazard and spatial unit data.
        """
        self.hazard_data = None
        self.spatial_units = None
        self.pop = None

    def prep_data(self, path_to_data: str, geo_type: str) -> gpd.GeoDataFrame:
        """
        Read, clean, and preprocess geospatial data for for PopEstimator.

        This method loads a geospatial file (GeoJSON or Parquet) of either hazard
        data or spatial unit data and prepares the data to be passed to
        estimate_exposed_pop. It cleans and validates geometries in these
        files, ensuring they are not null or empty, and reprojects them to a
        consistent CRS. For hazard data, it creates buffered geometries based
        on user-defined buffer distance columns. The processed data is stored
        as a class attribute and returned.

        For hazard data:
        - The input must have a string column 'ID_hazard', a geometry column
        'geometry', and one or more numeric columns starting with 'buffer_dist_'
        (buffer distances in meters), and a coordinate reference system (CRS),
        but a specific CRS is not required.
        - Geometries are cleaned (null/empty removed, made valid) and
        reprojected to WGS84.
        - Each geometry is assigned its best UTM projection (based on centroid).
        - For each buffer distance, a new column 'buffered_hazard_{suffix}' is
        created with the buffered geometry. Hazard geometries are first
        reprojected to the best UTM projection for that geometry, buffered in
        that projection, and then reprojected back to the original projection
        of the input data (WGS84).
        - The returned GeoDataFrame contains 'ID_hazard' and all buffered hazard
        columns.

        For spatial unit data:
        - The input must have a string column 'ID_spatial_unit' and a geometry
        column 'geometry', and a coordinate reference system (CRS), but a
        specific CRS is not required.
        - Geometries are cleaned (null/empty removed, made valid) and reprojected
        to WGS84.
        - The returned GeoDataFrame contains all columns from the input.

        This method sets PopEstimator class attributes 'hazards' and 'spatial_units'
        and returns the processed GeoDataFrame.
        If the file is empty, returns None and sets attributes to None.

        Parameters
        ----------
        path_to_data : str
            Path to the input geospatial data file (.geojson or .parquet) with
            required columns:
            'ID_hazard' or 'ID_spatial_unit', 'geometry', and buffer distance
            columns (for hazard data).
        geo_type : str
            Type of data to process: 'hazard' or 'spatial_unit'.

        Returns
        -------
        geopandas.GeoDataFrame or None
        Cleaned and processed GeoDataFrame with ID columns, geometries, and
        (for hazards) buffered hazard columns.
        """
        shp_df = self._read_data(path_to_data)
        if shp_df.empty:
            return None

        shp_df = self._remove_missing_geometries(shp_df)
        shp_df = self._make_geometries_valid(shp_df)
        shp_df = self._reproject_to_wgs84(shp_df)

        if geo_type == "hazard":
            shp_df = self._add_utm_projection(shp_df)
            shp_df = self._add_buffered_geoms(shp_df)
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

        elif geo_type == "spatial_unit":
            self.spatial_units = shp_df
            return shp_df

        else:
            raise ValueError("geo_type must be 'hazard' or 'spatial_unit'")

    def exposed_pop(
        self,
        pop_path: str,
        hazard_specific: bool,
        hazards: gpd.GeoDataFrame = None,
        spatial_units: gpd.GeoDataFrame = None,
    ) -> pd.DataFrame:
        """
        Estimate the population exposed to hazards, optionally within spatial units.

        This method calculates the sum of raster values (e.g., population) within hazard
        geometries, or within the intersection of hazard geometries and spatial units.
        It supports both hazard-specific and combined hazard analyses, and can use
        pre-loaded or provided hazard and spatial unit data.

        There are four main cases, depending on the combination of arguments:

        1. **Hazard-specific exposure, no spatial units** (`hazard_specific=True`, `spatial_units=None`):
            - Calculates the exposed population for each hazard geometry (and each buffer, if present).
            - Returns a DataFrame with one row per hazard and one or more 'exposed' columns.

        2. **Combined hazards, no spatial units** (`hazard_specific=False`, `spatial_units=None`):
            - All hazard geometries (for each buffer) are merged into a single geometry.
            - Calculates the total exposed population for the union of all hazards.
            - Returns a DataFrame with a single row and one or more 'exposed' columns.

        3. **Hazard-specific exposure within spatial units** (`hazard_specific=True`, `spatial_units` provided):
            - Calculates the exposed population for each intersection of hazard geometry (and buffer) with each spatial unit.
            - Returns a DataFrame with one row per (hazard, spatial unit) pair and one or more 'exposed' columns.

        4. **Combined hazards within spatial units** (`hazard_specific=False`, `spatial_units` provided):
            - All hazard geometries (for each buffer) are merged into a single geometry.
            - Calculates the exposed population for the intersection of the combined hazard geometry with each spatial unit.
            - Returns a DataFrame with one row per spatial unit and one or more 'exposed' columns.

        The method can use hazards and spatial units provided as arguments, or use those previously loaded
        via `prepare_data`. All geometries are cleaned and reprojected as needed. The population raster
        should be provided as a file path and must have a CRS. The method ensures CRS alignment between
        vector and raster data, and ignores invalid or empty geometries.

        The resulting DataFrame contains the relevant ID columns (`ID_hazard`, `ID_spatial_unit`) and
        one or more 'exposed' columns, corresponding to the sum of raster values within each geometry
        or intersection. Geometry columns are not included in the output.

        This method sets the class attribute `exposed` and returns the DataFrame.

        Parameters
        ----------
        pop_path : str
            Path to the population raster file.
        hazard_specific : bool
            If True, exposure is calculated for each hazard individually.
            If False, hazard geometries are combined before exposure calculation.
        hazards : geopandas.GeoDataFrame, optional
            GeoDataFrame containing hazard geometries and buffer columns.
            If None, uses self.hazard_data.
        spatial_units : geopandas.GeoDataFrame, optional
            GeoDataFrame containing spatial unit geometries.
            If None, uses self.spatial_units.

        Returns
        -------
        pandas.DataFrame
            DataFrame with ID columns and one or more 'exposed' columns indicating
            the sum of raster values within each geometry or intersection.
            The DataFrame does not include geometry columns.
        """

        if hazards is None:
            hazards = self.hazard_data
        if spatial_units is None:
            spatial_units = self.spatial_units
        if hazards is None:
            return None

        if spatial_units is None:
            if not hazard_specific:
                hazards = self._combine_geometries(hazards)
            exposed = self._mask_raster_partial_pixel(hazards, pop_path)
            self.exposed = exposed
            return exposed

        else:
            if not hazard_specific:
                hazards = self._combine_geometries(hazards)
            intersected_hazards = self._get_unit_hazard_intersections(
                hazards=hazards, spatial_units=spatial_units
            )
            exposed = self._mask_raster_partial_pixel(
                intersected_hazards, raster_path=pop_path
            )
            self.exposed = exposed
            return exposed

    def pop(self, pop_path: str, spatial_units: str) -> pd.DataFrame:
        """
        Estimate the total population residing within each spatial unit according
        to a population raster.

        This method calculates the sum of raster values (e.g., population)
        within the boundaries of each spatial unit geometry provided. It uses
        the exact_extract method to perform partial-pixel extraction from the
        raster, ensuring accurate population estimates even when spatial unit
        boundaries do not align perfectly with raster cells.

        The method expects a GeoDataFrame of spatial units, each with a unique
        'ID_spatial_unit' and a valid geometry column called 'geometry'. The
        population raster should be provided as a file path and have a CRS.
        The method will reproject the spatial units to match the CRS of the
        raster if necessary, and will ignore any spatial units with null or
        invalid geometries.

        The resulting DataFrame contains the spatial unit ID column and a
        population column called 'population'.

        This method is meant to be used with the same population raster as
        'estimate_exposed_pop' to provide denominators for the total population
        in each spatial unit.

        This method sets the class attribute pop and returns a dataframe.

        Parameters
        ----------
        pop_path : str
            Path to the population raster file (.tif).
        spatial_units : geopandas.GeoDataFrame
            GeoDataFrame containing spatial unit geometries. Must include a column
            'ID_spatial_unit' and a valid geometry column called 'geometry'.

        Returns
        -------
        pandas.DataFrame
            DataFrame with 'ID_spatial_unit' and a 'population' column,
            where each value is the sum of raster values within the
            corresponding spatial unit.
            The DataFrame does not include geometry columns.
        """
        residing = self._mask_raster_partial_pixel(spatial_units, raster_path=pop_path)
        residing = residing.rename(
            columns=lambda c: c.replace("exposedgeometry", "population")
        )
        self.pop = residing
        return residing

    # --- Helper methods below ---

    def _read_data(self, path: str) -> gpd.GeoDataFrame:
        """
        Read geospatial data for PopEstimator from a file into a GeoDataFrame.

        This method supports both .geojson and .parquet file formats, but hazard
        data must have a str column called 'ID_hazard', numeric columns starting
        with 'buffer_dist_' for buffer distances, and a geometry column called
        'geometry'. Spatial unit data must have a str column called
        'ID_spatial_unit' and a geometry column called 'geometry'.

        :param path: Path to the data file (.geojson or .parquet).
        :type path: str
        :returns: Loaded GeoDataFrame.
        :rtype: geopandas.GeoDataFrame
        :raises FileNotFoundError: If the file type is unsupported.
        """
        path = Path(path)
        if path.suffix == ".geojson":
            shp_df = gpd.read_file(path)
        elif path.suffix == ".parquet":
            shp_df = gpd.read_parquet(path)
        else:
            raise FileNotFoundError(f"File not found or unsupported file type: {path}")
        return shp_df

    def _remove_missing_geometries(self, shp_df: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Remove rows from hazard dataframe or spatial unit dataframe with null
        or empty geometries.

        :param shp_df: Input GeoDataFrame.
        :type shp_df: geopandas.GeoDataFrame
        :returns: GeoDataFrame with missing geometries removed.
        :rtype: geopandas.GeoDataFrame
        """
        return shp_df[shp_df["geometry"].notnull() & ~shp_df["geometry"].is_empty]

    def _make_geometries_valid(self, shp_df: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Make all geometries in the GeoDataFrame valid.

        :param shp_df: Input GeoDataFrame.
        :type shp_df: geopandas.GeoDataFrame
        :returns: GeoDataFrame with valid geometries.
        :rtype: geopandas.GeoDataFrame
        """
        shp_df["geometry"] = shp_df["geometry"].apply(make_valid)
        return shp_df

    def _reproject_to_wgs84(self, shp_df: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Reproject all geometries in the GeoDataFrame to WGS84 (EPSG:4326).

        :param shp_df: Input GeoDataFrame.
        :type shp_df: geopandas.GeoDataFrame
        :returns: GeoDataFrame with geometries reprojected to WGS84.
        :rtype: geopandas.GeoDataFrame
        """
        if shp_df.crs != "EPSG:4326":
            shp_df = shp_df.to_crs("EPSG:4326")
        return shp_df

    def _get_best_utm_projection(self, lat, lon):
        """
        Return a string representation of the UTM projection EPSG code for a
        given latitude and longitude.

        :param lat: Latitude.
        :type lat: float
        :param lon: Longitude.
        :type lon: float
        :returns: EPSG code string for the best UTM projection.
        :rtype: str
        """
        zone_number = (lon + 180) // 6 + 1
        hemisphere = 326 if lat >= 0 else 327
        epsg_code = hemisphere * 100 + zone_number
        return f"EPSG:{int(epsg_code)}"

    def _add_utm_projection(self, shp_df: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Add a column with the best UTM projection for each geometry to the
        inputted GeoDataFrame.

        :param ch_shp: Input GeoDataFrame.
        :type ch_shp: geopandas.GeoDataFrame
        :returns: GeoDataFrame with UTM projection column added.
        :rtype: geopandas.GeoDataFrame
        """
        # get geom lat and lon
        shp_df["centroid_lon"] = shp_df.centroid.x
        shp_df["centroid_lat"] = shp_df.centroid.y
        # get projection for each hazard
        shp_df["utm_projection"] = shp_df.apply(
            lambda row: self._get_best_utm_projection(
                lat=row["centroid_lat"], lon=row["centroid_lon"]
            ),
            axis=1,
        )
        # select id, geometry, buffer dist, and utm projection
        buffer_cols = [col for col in shp_df.columns if col.startswith("buffer_dist")]
        shp_df = shp_df[["ID_hazard"] + buffer_cols + ["geometry", "utm_projection"]]
        return shp_df

    def _get_buffered_geom(self, row, buffer_col):
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

    def _add_buffered_geoms(self, shp_df: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Add one column for each buffer_dist_ column in the inputted dataframe
        containing buffered geometries, buffered with each buffer distance passed.
        The buffered geometries are named 'buffered_hazard_' + buffer distance
        column name.

        For example, if the inputted dataframe has a column called
        'buffer_dist_100', the outputted dataframe will have a column called
        'buffered_hazard_100'.
        The buffer distance is in meters.
        The buffered geometries are created in the best UTM projection for each
        geometry, and then reprojected back to the original projection of the
        inputted dataframe.

        :param shp_df: Input GeoDataFrame.
        :type shp_df: geopandas.GeoDataFrame
        :returns: GeoDataFrame with buffered hazard geometry column added.
        :rtype: geopandas.GeoDataFrame
        """
        buffer_cols = [col for col in shp_df.columns if col.startswith("buffer_dist")]
        for buffer_col in buffer_cols:
            suffix = buffer_col.replace("buffer_dist", "").strip("_")
            new_col = f"buffered_hazard_{suffix}" if suffix else "buffered_hazard"
            shp_df[new_col] = shp_df.apply(
                lambda row: self._get_buffered_geom(row, buffer_col), axis=1
            )

        return shp_df

    def _combine_geometries(
        self,
        shp_df: gpd.GeoDataFrame,
    ) -> gpd.GeoDataFrame:
        """
        Combine all geometries in columns starting with 'buffered_hazard' into a
        single geometry. Use chunks for efficiency.

        :param shp_df: Input GeoDataFrame.
        :type shp_df: geopandas.GeoDataFrame
        :returns: GeoDataFrame with one row and merged geometry columns for
        each buffer containing merged geometries.
        :rtype: geopandas.GeoDataFrame
        """
        chunk_size = 500
        buffered_cols = [
            col for col in shp_df.columns if col.startswith("buffered_hazard")
        ]
        merged_geoms = {}
        for col in buffered_cols:
            geoms = [
                g
                for g in shp_df[col]
                if g is not None and g.is_valid and not g.is_empty
            ]
            chunks = [
                geoms[i : i + chunk_size] for i in range(0, len(geoms), chunk_size)
            ]
            partial_unions = [unary_union(chunk) for chunk in chunks]
            final_union = unary_union(partial_unions)
            merged_geoms[col] = [final_union]
        merged_geoms["ID_hazard"] = ["merged_geoms"]
        combined_gdf = gpd.GeoDataFrame(
            merged_geoms, geometry=buffered_cols[0], crs=shp_df.crs
        )
        return combined_gdf

    def _get_unit_hazard_intersections(self, hazards, spatial_units):
        intersections = {}
        for col in [c for c in hazards.columns if c.startswith("buffered_hazard")]:
            # Select only ID_hazard and the current geometry column
            hazards_subset = hazards[["ID_hazard", col]].copy()
            hazards_geom = hazards_subset.set_geometry(col, crs=hazards.crs)
            intersection = gpd.overlay(hazards_geom, spatial_units, how="intersection")
            intersection = self._remove_missing_geometries(intersection)
            intersection = self._make_geometries_valid(intersection)
            intersection = intersection.rename_geometry(col)
            intersection = intersection.set_geometry(col, crs=hazards.crs)

            intersections[col] = intersection
        intersected_dfs = [
            df for df in intersections.values() if df is not None and not df.empty
        ]

        intersected_hazards = functools.reduce(
            lambda left, right: pd.merge(
                left, right, on=["ID_hazard", "ID_spatial_unit"], how="outer"
            ),
            intersected_dfs,
        )
        return intersected_hazards

    def _mask_raster_partial_pixel(self, shp_df: gpd.GeoDataFrame, raster_path: str):
        """
        Calculate the sum of raster values (e.g., population) within each
        geometry using exact_extract.

        :param shp_df: Input GeoDataFrame.
        :type shp_df: geopandas.GeoDataFrame
        :param raster_path: Path to the raster file.
        :type raster_path: str
        :returns: GeoDataFrame with an 'exposed' column containing the sum for each geometry.
        :rtype: geopandas.GeoDataFrame
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
            if col.startswith("exposed") or col in ["ID_hazard", "ID_spatial_unit"]
        ]
        shp_exposed = shp_df[cols]

        return shp_exposed
