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

from popexposure.pop_estimator import PopEstimator


class TestPopEstimator:
    """Test cases for PopEstimator methods."""

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
    def temp_dir(self):
        """Create a temporary directory for test files."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)

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

    @pytest.fixture
    def pop_estimator(self, standard_raster, quadrant_admin_units, temp_dir):
        """Create a PopEstimator instance for testing with required pop and admin data."""
        # Save admin units to file (GeoJSON)
        admin_path = os.path.join(temp_dir, "admin.geojson")
        quadrant_admin_units.to_file(admin_path, driver="GeoJSON")
        # Pass raster path and admin path to PopEstimator
        return PopEstimator(pop_data=standard_raster, admin_data=admin_path)

    @pytest.fixture
    def pop_estimator_no_admin(self, standard_raster):
        """Create a PopEstimator instance for testing without admin data."""
        return PopEstimator(pop_data=standard_raster)

    def save_gdf(self, gdf, temp_dir, filename):
        """Helper to save GeoDataFrame to file."""
        path = os.path.join(temp_dir, filename)
        gdf.to_file(path, driver="GeoJSON")
        return path

    def test_prep_data_removes_missing_geometries(self, pop_estimator, temp_dir):
        """Test _process_hazard_data removes null geometries and creates buffered columns."""
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
        result = pop_estimator._process_hazard_data(path)

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
        self, pop_estimator_no_admin, temp_dir, point_hazards
    ):
        """Test exposure calculation with point hazards."""
        hazard_path = self.save_gdf(point_hazards, temp_dir, "hazards.geojson")

        result = pop_estimator_no_admin.est_exposed_pop(
            hazard_specific=True, hazard_data=hazard_path
        )
        print(result)

        assert len(result) == 3
        assert set(result.columns) == {"ID_hazard", "exposed_1", "exposed_10"}
        assert set(result["ID_hazard"]) == {"point_a", "point_b", "point_c"}

        # Check expected exposure values
        for _, row in result.iterrows():
            assert abs(row["exposed_1"] - 0.00032) < 0.01
            assert abs(row["exposed_10"] - 0.0032) < 0.01

    def test_est_exposed_pop_with_admin_units(self, pop_estimator, temp_dir):
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

        print(center_hazard)
        center_hazard = center_hazard.set_geometry("geometry")

        hazard_path = self.save_gdf(center_hazard, temp_dir, "hazard.geojson")

        result = pop_estimator.est_exposed_pop(
            hazard_specific=True,
            hazard_data=hazard_path,
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
        self, pop_estimator_no_admin, temp_dir, point_hazards
    ):
        """Test hazard_specific=False combines all hazards."""
        hazard_path = self.save_gdf(point_hazards, temp_dir, "hazards.geojson")

        result = pop_estimator_no_admin.est_exposed_pop(
            hazard_specific=False, hazard_data=hazard_path
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
        self, pop_estimator, temp_dir, quadrant_admin_units
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

        result = pop_estimator.est_exposed_pop(
            hazard_specific=False, hazard_data=hazard_path
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

        # For this test, we need a new PopEstimator with the new raster and admin units
        pop_estimator_area = PopEstimator(pop_data=raster_path, admin_data=admin_path)
        result = pop_estimator_area.est_total_pop()

        assert len(result) == 4
        assert set(result.columns) == {"ID_admin_unit", "population"}

        # Each 10x10 square should have population ≈ 100
        for _, row in result.iterrows():
            assert abs(row["population"] - 100.0) < 10

    def test_is_admin_data_prepped(self, pop_estimator, quadrant_admin_units, temp_dir):
        """Test _is_admin_data_prepped method with processed and unprocessed data."""
        # Test with unprocessed data (should return False)
        unprocessed_admin = quadrant_admin_units.copy()
        # Change CRS to make it unprocessed
        unprocessed_admin = unprocessed_admin.to_crs("EPSG:3857")
        assert not pop_estimator._is_admin_data_prepped(unprocessed_admin)

        # Test with data missing ID column (should return False)
        no_id_admin = quadrant_admin_units.copy()
        no_id_admin = no_id_admin.rename(columns={"ID_admin_unit": "name"})
        assert not pop_estimator._is_admin_data_prepped(no_id_admin)

        # Test with processed data (should return True)
        processed_admin = pop_estimator._process_admin_data(quadrant_admin_units)
        assert pop_estimator._is_admin_data_prepped(processed_admin)

        # Test with null geometries (should return False)
        null_geom_admin = quadrant_admin_units.copy()
        null_geom_admin.loc[0, "geometry"] = None
        assert not pop_estimator._is_admin_data_prepped(null_geom_admin)

    def test_is_hazard_data_prepped(self, pop_estimator, point_hazards, temp_dir):
        """Test _is_hazard_data_prepped method with processed and unprocessed data."""
        # Test with unprocessed data (should return False)
        unprocessed_hazards = point_hazards.copy()
        assert not pop_estimator._is_hazard_data_prepped(unprocessed_hazards)

        # Test with data missing ID_hazard column (should return False)
        no_id_hazards = point_hazards.copy()
        no_id_hazards = no_id_hazards.rename(columns={"ID_hazard": "name"})
        assert not pop_estimator._is_hazard_data_prepped(no_id_hazards)

        # Test with processed data (should return True)
        processed_hazards = pop_estimator._process_hazard_data(point_hazards)
        assert pop_estimator._is_hazard_data_prepped(processed_hazards)

        # Test with wrong CRS (should return False)
        wrong_crs_hazards = processed_hazards.copy()
        wrong_crs_hazards = wrong_crs_hazards.to_crs("EPSG:3857")
        assert not pop_estimator._is_hazard_data_prepped(wrong_crs_hazards)

        # Test with missing buffered columns (should return False)
        no_buffer_hazards = point_hazards.copy()
        no_buffer_hazards = no_buffer_hazards[["ID_hazard", "geometry"]]
        assert not pop_estimator._is_hazard_data_prepped(no_buffer_hazards)

        # Test with geometry not set to buffered column (should return False)
        wrong_geom_hazards = processed_hazards.copy()
        # Create a dummy geometry column with valid geometries but wrong name
        wrong_geom_hazards["dummy_geom"] = wrong_geom_hazards.geometry
        wrong_geom_hazards = wrong_geom_hazards.set_geometry("dummy_geom")
        assert not pop_estimator._is_hazard_data_prepped(wrong_geom_hazards)

    def test_est_exposed_pop_with_prepped_data(
        self, pop_estimator, point_hazards, temp_dir
    ):
        """Test that est_exposed_pop works with both prepped and unprepped hazard data."""
        # Test with unprocessed data (file path)
        hazard_path = self.save_gdf(point_hazards, temp_dir, "hazards.geojson")
        result_unprocessed = pop_estimator.est_exposed_pop(
            hazard_data=hazard_path, hazard_specific=True
        )

        # Test with prepped GeoDataFrame
        processed_hazards = pop_estimator._process_hazard_data(point_hazards)
        result_prepped = pop_estimator.est_exposed_pop(
            hazard_data=processed_hazards, hazard_specific=True
        )

        # Results should be identical
        assert len(result_unprocessed) == len(result_prepped)
        assert set(result_unprocessed.columns) == set(result_prepped.columns)

        # Sort both dataframes by ID_hazard for comparison
        result_unprocessed = result_unprocessed.sort_values("ID_hazard").reset_index(
            drop=True
        )
        result_prepped = result_prepped.sort_values("ID_hazard").reset_index(drop=True)

        # Compare values (allowing for small floating point differences)
        for col in result_unprocessed.columns:
            if col not in ["ID_hazard", "ID_admin_unit"]:  # Skip string columns
                assert np.allclose(
                    result_unprocessed[col], result_prepped[col], rtol=1e-10
                )
            else:
                assert result_unprocessed[col].equals(result_prepped[col])

    def test_constructor_with_prepped_admin_data(
        self, standard_raster, quadrant_admin_units, temp_dir
    ):
        """Test that constructor works with both prepped and unprepped admin data."""
        # Create estimator with file path (unprocessed)
        admin_path = self.save_gdf(quadrant_admin_units, temp_dir, "admin.geojson")
        estimator_file = PopEstimator(pop_data=standard_raster, admin_data=admin_path)

        # Create estimator with prepped GeoDataFrame
        processed_admin = PopEstimator(pop_data=standard_raster)._process_admin_data(
            quadrant_admin_units
        )
        estimator_prepped = PopEstimator(
            pop_data=standard_raster, admin_data=processed_admin
        )

        # Both should have admin_data with same content
        assert estimator_file.admin_data is not None
        assert estimator_prepped.admin_data is not None
        assert len(estimator_file.admin_data) == len(estimator_prepped.admin_data)
        assert set(estimator_file.admin_data.columns) == set(
            estimator_prepped.admin_data.columns
        )

        # Test that both work for exposure calculation
        test_hazard = gpd.GeoDataFrame(
            {
                "ID_hazard": ["test"],
                "buffer_dist_10": [10],
                "geometry": [Point(-121.5, 37.5)],
            },
            crs="EPSG:4326",
        )

        result_file = estimator_file.est_exposed_pop(
            hazard_data=test_hazard, hazard_specific=True
        )
        result_prepped = estimator_prepped.est_exposed_pop(
            hazard_data=test_hazard, hazard_specific=True
        )

        # Results should be identical
        assert len(result_file) == len(result_prepped)
        for col in result_file.columns:
            if col not in ["ID_hazard", "ID_admin_unit"]:  # Skip string columns
                assert np.allclose(result_file[col], result_prepped[col], rtol=1e-10)

    def test_est_exposed_pop_prepped_vs_unprepped_comprehensive(
        self, pop_estimator, point_hazards, temp_dir
    ):
        """Test est_exposed_pop with prepped vs unprepped data across different scenarios."""
        # Test scenario 1: hazard_specific=True
        hazard_path = self.save_gdf(point_hazards, temp_dir, "hazards1.geojson")
        processed_hazards = pop_estimator._process_hazard_data(point_hazards)

        result_unprocessed_specific = pop_estimator.est_exposed_pop(
            hazard_data=hazard_path, hazard_specific=True
        )
        result_prepped_specific = pop_estimator.est_exposed_pop(
            hazard_data=processed_hazards, hazard_specific=True
        )

        # Sort both by ID columns for comparison
        result_unprocessed_specific = result_unprocessed_specific.sort_values(
            ["ID_hazard", "ID_admin_unit"]
        ).reset_index(drop=True)
        result_prepped_specific = result_prepped_specific.sort_values(
            ["ID_hazard", "ID_admin_unit"]
        ).reset_index(drop=True)

        assert len(result_unprocessed_specific) == len(result_prepped_specific)
        assert set(result_unprocessed_specific.columns) == set(
            result_prepped_specific.columns
        )

        # Compare numeric columns
        for col in result_unprocessed_specific.columns:
            if col not in ["ID_hazard", "ID_admin_unit"]:
                assert np.allclose(
                    result_unprocessed_specific[col],
                    result_prepped_specific[col],
                    rtol=1e-10,
                ), f"Column {col} differs between prepped and unprepped data"

        # Test scenario 2: hazard_specific=False
        result_unprocessed_combined = pop_estimator.est_exposed_pop(
            hazard_data=hazard_path, hazard_specific=False
        )
        result_prepped_combined = pop_estimator.est_exposed_pop(
            hazard_data=processed_hazards, hazard_specific=False
        )

        # Sort both by ID columns for comparison
        result_unprocessed_combined = result_unprocessed_combined.sort_values(
            ["ID_admin_unit"]
        ).reset_index(drop=True)
        result_prepped_combined = result_prepped_combined.sort_values(
            ["ID_admin_unit"]
        ).reset_index(drop=True)

        assert len(result_unprocessed_combined) == len(result_prepped_combined)
        assert set(result_unprocessed_combined.columns) == set(
            result_prepped_combined.columns
        )

        # Compare numeric columns
        for col in result_unprocessed_combined.columns:
            if col not in ["ID_hazard", "ID_admin_unit"]:
                assert np.allclose(
                    result_unprocessed_combined[col],
                    result_prepped_combined[col],
                    rtol=1e-10,
                ), f"Column {col} differs between prepped and unprepped data in combined mode"

    def test_est_exposed_pop_prepped_vs_unprepped_no_admin(
        self, pop_estimator_no_admin, point_hazards, temp_dir
    ):
        """Test est_exposed_pop with prepped vs unprepped data without admin units."""
        hazard_path = self.save_gdf(point_hazards, temp_dir, "hazards2.geojson")
        processed_hazards = pop_estimator_no_admin._process_hazard_data(point_hazards)

        # Test scenario 1: hazard_specific=True, no admin
        result_unprocessed_specific = pop_estimator_no_admin.est_exposed_pop(
            hazard_data=hazard_path, hazard_specific=True
        )
        result_prepped_specific = pop_estimator_no_admin.est_exposed_pop(
            hazard_data=processed_hazards, hazard_specific=True
        )

        # Sort both by ID_hazard for comparison
        result_unprocessed_specific = result_unprocessed_specific.sort_values(
            "ID_hazard"
        ).reset_index(drop=True)
        result_prepped_specific = result_prepped_specific.sort_values(
            "ID_hazard"
        ).reset_index(drop=True)

        assert len(result_unprocessed_specific) == len(result_prepped_specific)
        assert set(result_unprocessed_specific.columns) == set(
            result_prepped_specific.columns
        )

        for col in result_unprocessed_specific.columns:
            if col != "ID_hazard":
                assert np.allclose(
                    result_unprocessed_specific[col],
                    result_prepped_specific[col],
                    rtol=1e-10,
                )

        # Test scenario 2: hazard_specific=False, no admin
        result_unprocessed_combined = pop_estimator_no_admin.est_exposed_pop(
            hazard_data=hazard_path, hazard_specific=False
        )
        result_prepped_combined = pop_estimator_no_admin.est_exposed_pop(
            hazard_data=processed_hazards, hazard_specific=False
        )

        assert len(result_unprocessed_combined) == len(result_prepped_combined) == 1
        assert set(result_unprocessed_combined.columns) == set(
            result_prepped_combined.columns
        )

        for col in result_unprocessed_combined.columns:
            if col != "ID_hazard":
                assert np.allclose(
                    result_unprocessed_combined[col],
                    result_prepped_combined[col],
                    rtol=1e-10,
                )

    def test_est_total_pop_prepped_vs_unprepped_admin(
        self, standard_raster, quadrant_admin_units, temp_dir
    ):
        """Test est_total_pop with prepped vs unprepped admin data."""
        # Create estimator with unprepped admin data (file path)
        admin_path = self.save_gdf(
            quadrant_admin_units, temp_dir, "admin_total.geojson"
        )
        estimator_unprepped = PopEstimator(
            pop_data=standard_raster, admin_data=admin_path
        )

        # Create estimator with prepped admin data
        processed_admin = PopEstimator(pop_data=standard_raster)._process_admin_data(
            quadrant_admin_units
        )
        estimator_prepped = PopEstimator(
            pop_data=standard_raster, admin_data=processed_admin
        )

        # Test with default stat="sum"
        result_unprepped_sum = estimator_unprepped.est_total_pop()
        result_prepped_sum = estimator_prepped.est_total_pop()

        # Sort both by ID_admin_unit for comparison
        result_unprepped_sum = result_unprepped_sum.sort_values(
            "ID_admin_unit"
        ).reset_index(drop=True)
        result_prepped_sum = result_prepped_sum.sort_values(
            "ID_admin_unit"
        ).reset_index(drop=True)

        assert len(result_unprepped_sum) == len(result_prepped_sum)
        assert set(result_unprepped_sum.columns) == set(result_prepped_sum.columns)
        assert set(result_unprepped_sum.columns) == {"ID_admin_unit", "population"}

        # Compare population values
        assert np.allclose(
            result_unprepped_sum["population"],
            result_prepped_sum["population"],
            rtol=1e-10,
        ), "Population values differ between prepped and unprepped admin data"

        # Verify ID_admin_unit values are the same
        assert result_unprepped_sum["ID_admin_unit"].equals(
            result_prepped_sum["ID_admin_unit"]
        ), "ID_admin_unit values differ between prepped and unprepped admin data"

        # Test with stat="mean"
        result_unprepped_mean = estimator_unprepped.est_total_pop(stat="mean")
        result_prepped_mean = estimator_prepped.est_total_pop(stat="mean")

        # Sort both by ID_admin_unit for comparison
        result_unprepped_mean = result_unprepped_mean.sort_values(
            "ID_admin_unit"
        ).reset_index(drop=True)
        result_prepped_mean = result_prepped_mean.sort_values(
            "ID_admin_unit"
        ).reset_index(drop=True)

        assert len(result_unprepped_mean) == len(result_prepped_mean)
        assert np.allclose(
            result_unprepped_mean["population"],
            result_prepped_mean["population"],
            rtol=1e-10,
        ), "Mean population values differ between prepped and unprepped admin data"

    def test_est_exposed_pop_different_stats_prepped_vs_unprepped(
        self, pop_estimator, point_hazards, temp_dir
    ):
        """Test est_exposed_pop with different stat parameters using prepped vs unprepped data."""
        hazard_path = self.save_gdf(point_hazards, temp_dir, "hazards_stats.geojson")
        processed_hazards = pop_estimator._process_hazard_data(point_hazards)

        # Test with stat="sum"
        result_unprepped_sum = pop_estimator.est_exposed_pop(
            hazard_data=hazard_path, hazard_specific=True, stat="sum"
        )
        result_prepped_sum = pop_estimator.est_exposed_pop(
            hazard_data=processed_hazards, hazard_specific=True, stat="sum"
        )

        # Sort both for comparison
        result_unprepped_sum = result_unprepped_sum.sort_values(
            ["ID_hazard", "ID_admin_unit"]
        ).reset_index(drop=True)
        result_prepped_sum = result_prepped_sum.sort_values(
            ["ID_hazard", "ID_admin_unit"]
        ).reset_index(drop=True)

        for col in result_unprepped_sum.columns:
            if col not in ["ID_hazard", "ID_admin_unit"]:
                assert np.allclose(
                    result_unprepped_sum[col], result_prepped_sum[col], rtol=1e-10
                ), f"Sum stat column {col} differs between prepped and unprepped data"

        # Test with stat="mean"
        result_unprepped_mean = pop_estimator.est_exposed_pop(
            hazard_data=hazard_path, hazard_specific=True, stat="mean"
        )
        result_prepped_mean = pop_estimator.est_exposed_pop(
            hazard_data=processed_hazards, hazard_specific=True, stat="mean"
        )

        # Sort both for comparison
        result_unprepped_mean = result_unprepped_mean.sort_values(
            ["ID_hazard", "ID_admin_unit"]
        ).reset_index(drop=True)
        result_prepped_mean = result_prepped_mean.sort_values(
            ["ID_hazard", "ID_admin_unit"]
        ).reset_index(drop=True)

        for col in result_unprepped_mean.columns:
            if col not in ["ID_hazard", "ID_admin_unit"]:
                assert np.allclose(
                    result_unprepped_mean[col], result_prepped_mean[col], rtol=1e-10
                ), f"Mean stat column {col} differs between prepped and unprepped data"
