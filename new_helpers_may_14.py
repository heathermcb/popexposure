import geopandas as gpd
import numpy as np
import pyarrow.parquet as pq
import rasterio
import rasterio.windows
from tqdm import tqdm
from shapely.validation import make_valid
from shapely.ops import unary_union
from shapely import wkt
from pathlib import Path
from exactextract import exact_extract
import warnings


class PopEstimator:

    def __init__(self):
        """
        Initialize the PopEstimator with empty attributes for hazard and spatial unit data.
        """
        self.hazard_data = None
        self.spatial_unit = None

    def prepare_data(self, path_to_data: str, geo_type: str) -> gpd.GeoDataFrame:
        """
        Read, clean, and preprocess geospatial data.

        :param path_to_data: Path to the input geospatial data file.
        :type path_to_data: str
        :param geo_type: Type of data to process ('hazard' or 'spatial_unit').
        :type geo_type: str
        :returns: Cleaned and processed GeoDataFrame.
        :rtype: geopandas.GeoDataFrame
        """
        self._print_geo_message(geo_type)
        shp_df = self._read_data(path_to_data)
        shp_df = self._remove_missing_geometries(shp_df)
        shp_df = self._make_geometries_valid(shp_df)
        shp_df = self._reproject_to_wgs84(shp_df)
        if geo_type == "hazard":
            shp_df = self._add_utm_projection(shp_df)
            shp_df = self._add_buffered_geom_col(shp_df)
            self.hazard_data = shp_df
        elif geo_type == "spatial_unit":
            self.spatial_unit = shp_df
        return shp_df

    def estimate_exposed_pop(
        self,
        pop_path: str,
        hazard_specific: bool,
        hazards: gpd.GeoDataFrame = None,
        spatial_unit: gpd.GeoDataFrame = None,
    ) -> gpd.GeoDataFrame:
        """
        Estimate the number of people exposed to hazards.

        :param pop_path: Path to the population raster file.
        :type pop_path: str
        :param hazard_specific: Whether to estimate exposure for each hazard separately.
        :type hazard_specific: bool
        :param hazards: GeoDataFrame of hazard geometries. If None, uses self.hazard_data.
        :type hazards: geopandas.GeoDataFrame, optional
        :param spatial_unit: GeoDataFrame of spatial units. If None, uses self.spatial_unit.
        :type spatial_unit: geopandas.GeoDataFrame, optional
        :returns: DataFrame with estimated exposed population.
        :rtype: geopandas.GeoDataFrame
        """
        if hazards is None:
            hazards = self.hazard_data
        if spatial_unit is None:
            spatial_unit = self.spatial_unit

        if hazards is None:
            raise ValueError("Hazard data must be provided or prepared first.")

        if spatial_unit is not None:
            hazards = hazards[["ID_hazard", "buffered_hazard"]]
            hazards = hazards.rename(columns={"buffered_hazard": "geometry"})
            hazards = hazards.set_geometry("geometry")
            if not hazard_specific:
                hazards = self._combine_geometries(hazards, id_column="ID_hazard")
            hazards = hazards.set_geometry("geometry", crs="EPSG:4326")
            unit_hazard_intersection = gpd.overlay(
                hazards, spatial_unit, how="intersection"
            )
            exposed = self._mask_raster_partial_pixel(
                unit_hazard_intersection, pop_path
            )
            exposed = exposed[["ID_hazard", "ID_spatial_unit", "exposed"]]
            self.exposed = exposed
            return exposed
        else:
            hazards = hazards[["ID_hazard", "buffered_hazard"]]
            hazards = hazards.rename(columns={"buffered_hazard": "geometry"})
            hazards = hazards.set_geometry("geometry", crs="EPSG:4326")
            if not hazard_specific:
                hazards = self._combine_geometries(hazards, id_column="ID_hazard")
            exposed = self._mask_raster_partial_pixel(hazards, pop_path)
            exposed = exposed[["ID_hazard", "exposed"]]
            self.exposed = exposed
            return exposed

        # ----------------------------------------------------------------------

        if not hazard_specific:
            if spatial_unit is not None:
                hazards = hazards[["ID_hazard", "buffered_hazard"]]
                hazards = hazards.rename(columns={"buffered_hazard": "geometry"})
                hazards = hazards.set_geometry("geometry")
                hazards = self._combine_geometries(hazards, id_column="ID_hazard")
                hazards = hazards.set_geometry("geometry", crs="EPSG:4326")
                unit_hazard_intersection = gpd.overlay(
                    hazards, spatial_unit, how="intersection"
                )
                exposed = self._mask_raster_partial_pixel(
                    unit_hazard_intersection, pop_path
                )
                exposed = exposed[["ID_hazard", "ID_spatial_unit", "exposed"]]
                self.exposed = exposed
                return exposed
            else:
                hazards = hazards[["ID_hazard", "buffered_hazard"]]
                hazards = hazards.rename(columns={"buffered_hazard": "geometry"})
                hazards = hazards.set_geometry("geometry", crs="EPSG:4326")
                hazards = self._combine_geometries(hazards, id_column="ID_hazard")
                exposed = self._mask_raster_partial_pixel(hazards, pop_path)
                exposed = exposed[["ID_hazard", "exposed"]]
                self.exposed = exposed
                return exposed

        else:
            if spatial_unit is not None:
                hazards = hazards[["ID_hazard", "buffered_hazard"]]
                hazards = hazards.rename(columns={"buffered_hazard": "geometry"})
                hazards = hazards.set_geometry("geometry")
                hazards = hazards.set_geometry("geometry", crs="EPSG:4326")
                unit_hazard_intersection = gpd.overlay(
                    hazards, spatial_unit, how="intersection"
                )
                exposed = self._mask_raster_partial_pixel(
                    unit_hazard_intersection, pop_path
                )
                exposed = exposed[["ID_hazard", "ID_spatial_unit", "exposed"]]
                self.exposed = exposed
                return exposed
            else:
                hazards = hazards[["ID_hazard", "buffered_hazard"]]
                hazards = hazards.rename(columns={"buffered_hazard": "geometry"})
                hazards = hazards.set_geometry("geometry", crs="EPSG:4326")
                exposed = self._mask_raster_partial_pixel(hazards, pop_path)
                exposed = exposed[["ID_hazard", "exposed"]]
                self.exposed = exposed
                return exposed

    # def estimate_pop(self, pop_path: str, spatial_unit: str):

    # --- Helper methods below ---

    def _print_geo_message(self, geo_type: str):
        """
        Print a message describing the type of geospatial data being processed.

        :param geo_type: Type of data ('hazard' or 'spatial_unit').
        :type geo_type: str
        """
        if geo_type == "hazard":
            print("Reading data and finding best UTM projection for hazard geometries")
        elif geo_type == "spatial_unit":
            print("Reading spatial unit geometries")

    def _read_data(self, path: str) -> gpd.GeoDataFrame:
        """
        Read geospatial data from a file.

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
        Remove rows from hazard dataframe or spatial unit dataframe with missing or empty geometries.

        :param shp_df: Input GeoDataFrame.
        :type shp_df: geopandas.GeoDataFrame
        :returns: GeoDataFrame with missing geometries removed.
        :rtype: geopandas.GeoDataFrame
        """
        return shp_df[~shp_df["geometry"].is_empty]

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
        if shp_df.crs != "EPSG:4326":
            shp_df = shp_df.to_crs("EPSG:4326")
        return shp_df

    def _get_best_utm_projection(self, lat, lon):
        """
        Calculate the best UTM projection EPSG code for a given latitude and longitude.

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
        Add a column with the best UTM projection for each geometry in the GeoDataFrame.

        :param ch_shp: Input GeoDataFrame.
        :type ch_shp: geopandas.GeoDataFrame
        :returns: GeoDataFrame with UTM projection column added.
        :rtype: geopandas.GeoDataFrame
        """
        # get lat and lon
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
        shp_df = shp_df[["ID_hazard", "buffer_dist", "geometry", "utm_projection"]]

        return shp_df

    def _add_buffered_geom_col(self, shp_df: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Add a column with buffered geometries for each row, using the best UTM projection.

        :param shp_df: Input GeoDataFrame.
        :type shp_df: geopandas.GeoDataFrame
        :returns: GeoDataFrame with buffered hazard geometry column added.
        :rtype: geopandas.GeoDataFrame
        """
        for index, row in tqdm(
            shp_df.iterrows(),
            total=len(shp_df),
            desc="Buffering hazard geometries",
        ):
            best_utm = row["utm_projection"]
            hazard_geom = row["geometry"]

            # create geoseries in best projection
            geom_series = gpd.GeoSeries([hazard_geom], crs=shp_df.crs)
            geom_series_utm = geom_series.to_crs(best_utm)

            # buffer distance is in meters
            buffer_dist = row["buffer_dist"]
            buffered_hazard_geometry = geom_series_utm.buffer(buffer_dist).iloc[0]
            # back to OG
            buffered_hazard_geometry = (
                gpd.GeoSeries([buffered_hazard_geometry], crs=best_utm)
                .to_crs(shp_df.crs)
                .iloc[0]
            )
            # add
            shp_df.at[index, "buffered_hazard"] = buffered_hazard_geometry

        return shp_df

    def _combine_geometries(
        self, shp_df: gpd.GeoDataFrame, id_column: str
    ) -> gpd.GeoDataFrame:
        """
        Combine all geometries in the GeoDataFrame into a single geometry, in chunks for efficiency.

        :param shp_df: Input GeoDataFrame.
        :type shp_df: geopandas.GeoDataFrame
        :param id_column: Name of the ID column to use for the merged geometry.
        :type id_column: str
        :returns: GeoDataFrame with a single merged geometry.
        :rtype: geopandas.GeoDataFrame
        """
        chunk_size = 500
        geoms = list(shp_df.geometry)
        chunks = [geoms[i : i + chunk_size] for i in range(0, len(geoms), chunk_size)]
        partial_unions = [unary_union(chunk) for chunk in chunks]
        final_union = unary_union(partial_unions)

        merged_id = "merged_geoms"
        combined_gdf = gpd.GeoDataFrame(
            {id_column: [merged_id], "geometry": [final_union]}, crs=shp_df.crs
        )
        return combined_gdf

    def _mask_raster_partial_pixel(self, shp_df: gpd.GeoDataFrame, raster_path: str):
        """
        Calculate the sum of raster values (e.g., population) within each geometry using exact_extract.

        :param shp_df: Input GeoDataFrame.
        :type shp_df: geopandas.GeoDataFrame
        :param raster_path: Path to the raster file.
        :type raster_path: str
        :returns: GeoDataFrame with an 'exposed' column containing the sum for each geometry.
        :rtype: geopandas.GeoDataFrame
        """
        print("Finding exposed population")
        # Open the raster file
        with rasterio.open(raster_path) as src:
            # Ensure CRS alignment
            if shp_df.crs != src.crs:
                shp_df = shp_df.to_crs(src.crs)

        # Use exact_extract to calculate population sums for each geometry
        # it returns a dictionary so we need to get out the right stuff
        num_people_affected = exact_extract(
            raster_path,
            shp_df,
            "sum",
        )
        sums = [hazard["properties"]["sum"] for hazard in num_people_affected]
        shp_df["exposed"] = sums
        # project shp_df back to wgs84
        shp_df = shp_df.to_crs("EPSG:4326")

        # final df
        return shp_df
