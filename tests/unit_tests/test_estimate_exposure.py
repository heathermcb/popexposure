"""
Unit tests for PopEstimator methods.
"""

import pytest
import geopandas as gpd
import pandas as pd
import tempfile
import os
import numpy as np
import rasterio
from rasterio.transform import from_bounds
from shapely.geometry import (
    Point,
    Polygon,
    LineString,
    MultiPolygon,
    GeometryCollection,
)
import sys
import shutil


# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from popexposure.estimate_exposure import PopEstimator


class TestPopEstimator:
    """Test cases for PopEstimator methods."""

    @pytest.fixture
    def pop_estimator(self):
        """Create a PopEstimator instance for testing."""
        return PopEstimator()

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)

    @pytest.fixture
    def standard_raster(self, temp_dir):
        """Create standard 100x100 raster with all 1s covering -122 to -121, 37 to 38."""
        raster_path = os.path.join(temp_dir, "raster_ones.tif")
        width, height = 100, 100
        west, south, east, north = -122.0, 37.0, -121.0, 38.0
        transform = from_bounds(west, south, east, north, width, height)
        raster_data = np.ones((height, width), dtype=np.float32)

        with rasterio.open(
            raster_path,
            "w",
            driver="GTiff",
            height=height,
            width=width,
            count=1,
            dtype=np.float32,
            crs="EPSG:4326",
            transform=transform,
            nodata=0,
        ) as dst:
            dst.write(raster_data, 1)
        return raster_path

    @pytest.fixture
    def point_hazards(self):
        """Create 3 point hazards with 1m and 10m buffers."""
        return gpd.GeoDataFrame(
            {
                "ID_hazard": ["point_a", "point_b", "point_c"],
                "buffer_dist_1": [1, 1, 1],
                "buffer_dist_10": [10, 10, 10],
                "geometry": [
                    Point(-121.5, 37.5),  # Center
                    Point(-121.3, 37.3),  # Southeast
                    Point(-121.7, 37.7),  # Northwest
                ],
            },
            crs="EPSG:4326",
        )

    @pytest.fixture
    def quadrant_admin_units(self):
        """Create 4 admin units in 2x2 grid meeting at center."""
        west, south, east, north = -122.0, 37.0, -121.0, 38.0
        center_lon, center_lat = -121.5, 37.5

        return gpd.GeoDataFrame(
            {
                "ID_admin_unit": ["northwest", "northeast", "southwest", "southeast"],
                "geometry": [
                    Polygon(
                        [
                            (west, center_lat),
                            (center_lon, center_lat),
                            (center_lon, north),
                            (west, north),
                            (west, center_lat),
                        ]
                    ),
                    Polygon(
                        [
                            (center_lon, center_lat),
                            (east, center_lat),
                            (east, north),
                            (center_lon, north),
                            (center_lon, center_lat),
                        ]
                    ),
                    Polygon(
                        [
                            (west, south),
                            (center_lon, south),
                            (center_lon, center_lat),
                            (west, center_lat),
                            (west, south),
                        ]
                    ),
                    Polygon(
                        [
                            (center_lon, south),
                            (east, south),
                            (east, center_lat),
                            (center_lon, center_lat),
                            (center_lon, south),
                        ]
                    ),
                ],
            },
            crs="EPSG:4326",
        )

    def save_gdf(self, gdf, temp_dir, filename):
        """Helper to save GeoDataFrame to file."""
        path = os.path.join(temp_dir, filename)
        gdf.to_file(path, driver="GeoJSON")
        return path

    def test_prep_data_removes_missing_geometries(self, pop_estimator, temp_dir):
        """Test prep_data removes null geometries and creates buffered columns."""
        # Create test data with one missing geometry
        hazards = gpd.GeoDataFrame(
            {
                "ID_hazard": ["point_1", "missing_1"],
                "buffer_dist_250": [250, 250],
                "buffer_dist_750": [750, 750],
                "geometry": [Point(-121.5, 37.5), None],
            },
            crs="EPSG:4326",
        )

        path = self.save_gdf(hazards, temp_dir, "hazards.geojson")
        result = pop_estimator.prep_data(path, "hazard")

        # Should have 1 row (missing geometry removed)
        assert len(result) == 1
        # Should have ID and buffered columns only
        assert set(result.columns) == {
            "ID_hazard",
            "buffered_hazard_250",
            "buffered_hazard_750",
        }
        # Should preserve the valid ID
        assert result["ID_hazard"].iloc[0] == "point_1"

    def test_est_exposed_pop_point_hazards(
        self, pop_estimator, temp_dir, point_hazards, standard_raster
    ):
        """Test exposure calculation with point hazards."""
        hazard_path = self.save_gdf(point_hazards, temp_dir, "hazards.geojson")
        prepared_hazards = pop_estimator.prep_data(hazard_path, "hazard")

        result = pop_estimator.est_exposed_pop(
            pop_path=standard_raster, hazard_specific=True, hazards=prepared_hazards
        )

        assert len(result) == 3
        assert set(result.columns) == {"ID_hazard", "exposed_1", "exposed_10"}
        assert set(result["ID_hazard"]) == {"point_a", "point_b", "point_c"}

        # Check expected exposure values
        for _, row in result.iterrows():
            assert abs(row["exposed_1"] - 0.00032) < 0.01
            assert abs(row["exposed_10"] - 0.0032) < 0.01

    def test_est_exposed_pop_with_admin_units(
        self, pop_estimator, temp_dir, quadrant_admin_units, standard_raster
    ):
        """Test exposure calculation with admin units using center point."""
        # Single center point hazard
        center_hazard = gpd.GeoDataFrame(
            {
                "ID_hazard": ["center_point"],
                "buffer_dist_10": [10],
                "geometry": [Point(-121.5, 37.5)],
            },
            crs="EPSG:4326",
        )

        hazard_path = self.save_gdf(center_hazard, temp_dir, "hazard.geojson")
        admin_path = self.save_gdf(quadrant_admin_units, temp_dir, "admin.geojson")

        prepared_hazards = pop_estimator.prep_data(hazard_path, "hazard")
        prepared_admin = pop_estimator.prep_data(admin_path, "admin_unit")

        result = pop_estimator.est_exposed_pop(
            pop_path=standard_raster,
            hazard_specific=True,
            hazards=prepared_hazards,
            admin_units=prepared_admin,
        )

        assert len(result) == 4
        assert set(result.columns) == {"ID_admin_unit", "ID_hazard", "exposed_10"}
        assert set(result["ID_admin_unit"]) == {
            "northwest",
            "northeast",
            "southwest",
            "southeast",
        }

        # Each admin unit should get exactly 1/4 of the total exposure
        # Total exposure for 10m buffer = 0.0032, so each quadrant = 0.0032/4 = 0.0008
        expected_exposure_per_quadrant = 0.0032 / 4  # = 0.0008

        for _, row in result.iterrows():
            assert (
                abs(row["exposed_10"] - expected_exposure_per_quadrant) < 0.01
            ), f"Admin unit {row['ID_admin_unit']} exposure {row['exposed_10']} should be ~{expected_exposure_per_quadrant}"

    def test_est_exposed_pop_combined_hazards(
        self, pop_estimator, temp_dir, point_hazards, standard_raster
    ):
        """Test hazard_specific=False combines all hazards."""
        hazard_path = self.save_gdf(point_hazards, temp_dir, "hazards.geojson")
        prepared_hazards = pop_estimator.prep_data(hazard_path, "hazard")

        result = pop_estimator.est_exposed_pop(
            pop_path=standard_raster, hazard_specific=False, hazards=prepared_hazards
        )

        assert len(result) == 1  # Combined into single row
        assert set(result.columns) == {
            "ID_hazard",
            "exposed_1",
            "exposed_10",
        }

        # Should be ~3x individual exposure (3 non-overlapping points)
        row = result.iloc[0]
        assert abs(row["exposed_1"] - 3 * 0.00032) < 0.01
        assert abs(row["exposed_10"] - 3 * 0.0032) < 0.01

    def test_est_exposed_pop_combined_with_admin_units(
        self, pop_estimator, temp_dir, quadrant_admin_units, standard_raster
    ):
        """Test hazard_specific=False with admin_units."""
        # Create 4 non-overlapping hazards, one per admin unit
        hazards = gpd.GeoDataFrame(
            {
                "ID_hazard": ["hazard_nw", "hazard_ne", "hazard_sw", "hazard_se"],
                "buffer_dist_10": [10, 10, 10, 10],
                "geometry": [
                    Point(-121.75, 37.75),  # NW center
                    Point(-121.25, 37.75),  # NE center
                    Point(-121.75, 37.25),  # SW center
                    Point(-121.25, 37.25),  # SE center
                ],
            },
            crs="EPSG:4326",
        )

        hazard_path = self.save_gdf(hazards, temp_dir, "hazards.geojson")
        admin_path = self.save_gdf(quadrant_admin_units, temp_dir, "admin.geojson")

        prepared_hazards = pop_estimator.prep_data(hazard_path, "hazard")
        prepared_admin = pop_estimator.prep_data(admin_path, "admin_unit")

        result = pop_estimator.est_exposed_pop(
            pop_path=standard_raster,
            hazard_specific=False,
            hazards=prepared_hazards,
            admin_units=prepared_admin,
        )
        assert len(result) == 4  # One row per admin unit
        assert set(result.columns) == {
            "ID_admin_unit",
            "ID_hazard",
            "exposed_10",
        }  # No ID_hazard
        assert set(result["ID_admin_unit"]) == {
            "northwest",
            "northeast",
            "southwest",
            "southeast",
        }

        # Each admin unit should have exactly 0.0032 exposure (one 10m point hazard each)
        expected_exposure = 0.0032

        for _, row in result.iterrows():
            assert (
                abs(row["exposed_10"] - expected_exposure) < 0.01
            ), f"Admin unit {row['ID_admin_unit']} exposure {row['exposed_10']} should be ~{expected_exposure}"

    def test_est_pop_returns_area(self, pop_estimator, temp_dir):
        """Test est_pop returns area when raster has all 1s."""
        # Create 4 simple 10x10 squares
        admin_units = gpd.GeoDataFrame(
            {
                "ID_admin_unit": ["admin_0_0", "admin_0_1", "admin_1_0", "admin_1_1"],
                "geometry": [
                    Polygon([(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]),
                    Polygon([(0, 10), (10, 10), (10, 20), (0, 20), (0, 10)]),
                    Polygon([(10, 0), (20, 0), (20, 10), (10, 10), (10, 0)]),
                    Polygon([(10, 10), (20, 10), (20, 20), (10, 20), (10, 10)]),
                ],
            },
            crs="EPSG:3857",
        )

        # Create matching raster (20x20 with all 1s)
        raster_path = os.path.join(temp_dir, "area_raster.tif")
        raster_data = np.ones((20, 20), dtype=np.float32)
        transform = rasterio.transform.from_bounds(0, 0, 20, 20, 20, 20)

        with rasterio.open(
            raster_path,
            "w",
            driver="GTiff",
            height=20,
            width=20,
            count=1,
            dtype=np.float32,
            crs="EPSG:3857",
            transform=transform,
        ) as dst:
            dst.write(raster_data, 1)

        admin_path = self.save_gdf(admin_units, temp_dir, "admin.geojson")
        prepared_admin = pop_estimator.prep_data(admin_path, "admin_unit")

        result = pop_estimator.est_pop(pop_path=raster_path, admin_units=prepared_admin)

        assert len(result) == 4
        assert set(result.columns) == {"ID_admin_unit", "population"}

        # Each 10x10 square should have population ≈ 100
        for _, row in result.iterrows():
            assert abs(row["population"] - 100.0) < 10
