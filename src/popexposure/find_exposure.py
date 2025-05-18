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
        Reads, cleans, and preprocesses geospatial data for exposure analysis.

        This function loads a geospatial file (GeoJSON or GeoParquet) containing either hazard data (e.g., wildfire burn zones, oil wells)
        or additional administrative geographies (e.g., ZCTAs, census tracts, referred to here as spatial units). It makes all geometries valid,
        removes empty geometries, and, for hazard data, generates buffered geometries for one or more user-specified buffer distances.
        Buffering is performed in the best Universal Transverse Mercator (UTM) projection based on each geometry's centroid latitude and longitude.

        Parameters
        ----------
        geo_type : str
            A string indicating the type of data to process. Must be either ``"hazard"`` for environmental hazard data or ``"spatial_unit"`` for administrative geography data.
        path_to_data : str
            Path to a geospatial data file (.geojson or .parquet). The file must contain either hazard data or administrative geography data, as specified by ``geo_type``.
            Data must have any coordinate reference system.
            - Hazard data must contain a string column ``"ID_hazard"`` with unique hazard IDs, a geometry column ``"geometry"``, and one or more numeric columns starting with ``"buffer_dist"`` with unique suffixes (e.g., ``"buffer_dist_main"``, ``"buffer_dist_1000"``) specifying buffer distances in meters. Buffer distances may be 0 or different for each hazard.
            - For spatial unit data, the file must contain a string column ``"ID_spatial_unit"`` with unique spatial unit IDs and a geometry column ``"geometry"``.

        Returns
        -------
        geopandas.GeoDataFrame or None
            A GeoDataFrame with cleaned and valid geometries.
            - If hazard data was passed, the output contains a column ``"ID_hazard"`` matching the input data, and one or more ``"buffered_hazard"`` geometry columns, with suffixes matching the passed ``buffer_dist`` columns (e.g., ``"buffered_hazard_main"``, ``"buffered_hazard_1000"``).
            - If spatial unit data was passed, the output contains columns ``"ID_spatial_unit"`` matching the input data and ``"geometry"``.
            - Empty geometries are removed.
            - If the input file is empty or contains no valid geometries, the function returns None.
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
        Estimate the number of people living within a buffer distance of environmental hazard(s) using a gridded population raster.

        This function calculates the sum of raster values within buffered hazard geometries, or within the intersection of buffered hazard geometries and additional administrative geographies, to find the population exposed to hazards. Users can choose to estimate either (a) hazard-specific counts (the number of people exposed to each unique buffered hazard in the set), or (b) a cumulative count (the number of unique people exposed to any of the input buffered hazards, avoiding double counting). Either estimate can be broken down by additional administrative geographies such as ZCTAs. Users must supply at least one buffered hazard column, but may supply additional buffered hazard columns to create estimates of exposure for different buffer distances.

        Parameters
        ----------
        hazards : geopandas.GeoDataFrame
            A GeoDataFrame with a coordinate reference system containing a string column called ``ID_hazard`` with unique hazard IDs, and one or more geometry columns starting with ``buffered_hazard`` containing buffered hazard geometries. ``buffered_hazard`` columns must each have a unique suffix (e.g., ``buffered_hazard_10``, ``buffered_hazard_100``, ``buffered_hazard_1000``).
        pop_path : str
            Path to a gridded population raster file, in TIFF format. The raster must have any coordinate reference system.
        hazard_specific : bool
            If True, exposure is calculated for each hazard individually (hazard-specific estimates). If False, geometries are combined before exposure is calculated, producing a single cumulative estimate.
        spatial_units : geopandas.GeoDataFrame, optional
            An optional GeoDataFrame of additional administrative geographies, containing a string column called ``ID_spatial_unit`` and a geometry column called ``geometry``.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with the following columns:
            - ``ID_hazard``: Always included.
            - ``ID_spatial_unit``: Included only if spatial units were provided.
            - One or more ``exposed`` columns: Each corresponds to a buffered hazard column (e.g., if the input had columns ``buffered_hazard_10``, ``buffered_hazard_100``, and ``buffered_hazard_1000``, the output will have ``exposed_10``, ``exposed_100``, and ``exposed_1000``). Each ``exposed`` column contains the sum of raster values (population) within the relevant buffered hazard geometry or buffered hazard geometry and spatial unit intersection.

            The number of rows in the output DataFrame depends on the function arguments:
            - If ``hazard_specific`` is True, the DataFrame contains one row per hazard or per hazard-spatial unit pair, if spatial units are provided.
            - If ``hazard_specific`` is False, the DataFrame contains a single row or one row per spatial unit, if spatial units are provided, with each ``exposed`` column representing the total population in the union of all buffered hazard geometries in that buffered hazard column.

        Notes
        -----
        There are four ways to use this function:

        1. Hazard-specific exposure, no additional administrative geographies (``hazard_specific=True, spatial_units=None``):
           Calculates the exposed population for each buffered hazard geometry. Returns a DataFrame with one row per hazard and one ``exposed`` column per buffered hazard column. If people lived within the buffer distance of more than one hazard, they are included in the exposure counts for each hazard they are near.

        2. Combined hazards, no additional administrative geographies (``hazard_specific=False, spatial_units=None``):
           All buffered hazard geometries in each buffered hazard column are merged into a single geometry, and the function calculates the total exposed population for the union of those buffered hazards. Returns a DataFrame with a single row and one ``exposed`` column for each buffered hazard column. If people were close to more than one hazard in the hazard set, they are counted once.

        3. Hazard-specific exposure within spatial units (``hazard_specific=True, spatial_units`` provided):
           Calculates the exposed population for each intersection of each buffered hazard geometry and each spatial unit. Returns a DataFrame with one row per buffered hazard-spatial unit pair and one ``exposed`` column per buffered hazard column. If people lived within the buffer distance of more than one hazard, they are included in the exposure counts for their spatial unit-hazard combination for each hazard they are near.

        4. Combined hazards within spatial units (``hazard_specific=False, spatial_units`` provided):
           All buffered hazard geometries in the same column are merged into a single geometry. Calculates the exposed population for the intersection of each buffered hazard combined geometry with each spatial unit. Returns a DataFrame with one row per spatial unit and one ``exposed`` column per buffered hazard column. If people were close to more than one hazard in the hazard set, they are counted once.
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
        Estimate the total population residing within administrative geographies using a gridded population raster.

        This function estimates the total population residing within administrative geographies (e.g., ZCTAs, census tracts) according to a provided gridded population raster. This method is meant to be used with the same population raster as ``exposed_pop`` to provide denominators for the total population in each administrative geography, allowing the user to compute the percentage of people exposed to hazards in each spatial unit. ``pop`` calculates the sum of raster values within the boundaries of each administrative geography geometry provided.

        Parameters
        ----------
        pop_path : str
            Path to a gridded population raster file, in TIFF format. The raster must have any coordinate reference system.
        spatial_units : geopandas.GeoDataFrame
            GeoDataFrame containing administrative geography geometries. Must include a string column called ``ID_spatial_unit`` with unique spatial unit IDs and a geometry column called ``geometry``.

        Returns
        -------
        pandas.DataFrame
            DataFrame with an ``ID_spatial_unit`` column matching the input and a ``population`` column, where each value is the sum of raster values within the corresponding spatial unit geometry.
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
