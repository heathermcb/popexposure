"""
Unit tests for reader.py

Creates test data with good and perverse cases, saves as different file formats,
and tests reading and validation functionality.
"""

import pytest
import geopandas as gpd
from pathlib import Path
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
import tempfile

import sys
from pathlib import Path

# Add the directory containing the module to Python path
module_path = Path("/Volumes/squirrel-utopia/popexposure/popexposure/utils")
sys.path.insert(0, str(module_path))

# Import the module directly
#import reader
from reader import *


class TestDataReader:
    """Test suite for DataReader class."""

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for test files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            yield Path(temp_dir)

    @pytest.fixture
    def hazard_test_data(self):
        """Create hazard data with good and perverse cases."""
        geometries = [
            # Good cases
            Point(-8235000, 4970000),  # Valid point
            Polygon(
                [
                    (-8237000, 4967000),
                    (-8235000, 4967000),
                    (-8235000, 4969000),
                    (-8237000, 4969000),
                    (-8237000, 4967000),
                ]
            ),  # Valid polygon
            # Perverse cases
            None,  # Empty geometry
            Polygon(
                [
                    (-8234000, 4967000),
                    (-8232000, 4969000),
                    (-8232000, 4967000),
                    (-8234000, 4969000),
                    (-8234000, 4967000),
                ]
            ),  # Self-intersecting bowtie
            LineString([(-8236000, 4968000), (-8234000, 4969000)]),  # Line geometry
        ]

        data = {
            "ID_hazard": ["h1", "h2", "h3", "h4", "h5"],
            "buffer_dist_500": [500, 750, 0, 1200, 800],
            "buffer_dist_1000": [1000, 1500, 0, 2400, 1600],
            "hazard_type": ["fire", "flood", "empty", "invalid", "line"],
            "geometry": geometries,
        }
        return gpd.GeoDataFrame(data, crs="EPSG:3857")

    @pytest.fixture
    def admin_unit_test_data(self):
        """Create admin unit data with good and perverse cases."""
        geometries = [
            # Good cases
            Polygon(
                [
                    (-8240000, 4970000),
                    (-8235000, 4970000),
                    (-8235000, 4975000),
                    (-8240000, 4975000),
                    (-8240000, 4970000),
                ]
            ),  # Valid polygon
            MultiPolygon(
                [
                    Polygon(
                        [
                            (-8245000, 4960000),
                            (-8240000, 4960000),
                            (-8240000, 4965000),
                            (-8245000, 4965000),
                            (-8245000, 4960000),
                        ]
                    ),
                    Polygon(
                        [
                            (-8243000, 4958000),
                            (-8242000, 4958000),
                            (-8242000, 4959000),
                            (-8243000, 4959000),
                            (-8243000, 4958000),
                        ]
                    ),
                ]
            ),  # Valid multipolygon
            # Perverse cases
            Point(-8237000, 4967000),  # Point instead of polygon
            LineString(
                [(-8238000, 4968000), (-8236000, 4968000)]
            ),  # Line instead of polygon
            None,  # Empty geometry
            Polygon(
                [
                    (-8234000, 4962000),
                    (-8232000, 4964000),
                    (-8232000, 4962000),
                    (-8234000, 4964000),
                    (-8234000, 4962000),
                ]
            ),  # Invalid self-intersecting
        ]

        data = {
            "ID_admin_unit": [
                "manhattan",
                "staten_island",
                "bad_point",
                "bad_line",
                "empty_unit",
                "invalid_polygon",
            ],
            "name": [
                "Manhattan",
                "Staten Island",
                "Bad Point",
                "Bad Line",
                "Empty Unit",
                "Invalid Polygon",
            ],
            "geometry": geometries,
        }
        return gpd.GeoDataFrame(data, crs="EPSG:3857")

    @pytest.fixture
    def test_files(self, hazard_test_data, admin_unit_test_data, temp_dir):
        """Create all test files: geojson, parquet, txt, and blank."""
        files = {}

        # Hazard files
        files["hazard_geojson"] = temp_dir / "hazards.geojson"
        files["hazard_parquet"] = temp_dir / "hazards.parquet"
        hazard_test_data.to_file(files["hazard_geojson"], driver="GeoJSON")
        hazard_test_data.to_parquet(files["hazard_parquet"])

        # Admin unit files
        files["admin_geojson"] = temp_dir / "admin_units.geojson"
        files["admin_parquet"] = temp_dir / "admin_units.parquet"
        admin_unit_test_data.to_file(files["admin_geojson"], driver="GeoJSON")
        admin_unit_test_data.to_parquet(files["admin_parquet"])

        # Unsupported files
        files["txt_file"] = temp_dir / "data.txt"
        files["txt_file"].write_text("This is not geospatial data")

        files["blank_file"] = temp_dir / "blank.geojson"
        files["blank_file"].touch()  # Create empty file

        return files

    # ===========================================
    # FILE READING TESTS
    # ===========================================

    def test_read_hazard_geojson(self, test_files):
        """Test reading hazard GeoJSON file."""
        gdf = read_geospatial_file(str(test_files["hazard_geojson"]))

        assert isinstance(gdf, gpd.GeoDataFrame)
        assert len(gdf) == 5
        assert "ID_hazard" in gdf.columns
        assert "geometry" in gdf.columns
        assert gdf.crs == "EPSG:3857"

    def test_read_hazard_parquet(self, test_files):
        """Test reading hazard Parquet file."""
        gdf = read_geospatial_file(str(test_files["hazard_parquet"]))

        assert isinstance(gdf, gpd.GeoDataFrame)
        assert len(gdf) == 5
        assert "ID_hazard" in gdf.columns
        assert "geometry" in gdf.columns
        assert gdf.crs == "EPSG:3857"

    def test_read_admin_geojson(self, test_files):
        """Test reading admin unit GeoJSON file."""
        gdf = read_geospatial_file(str(test_files["admin_geojson"]))

        assert isinstance(gdf, gpd.GeoDataFrame)
        assert len(gdf) == 6
        assert "ID_admin_unit" in gdf.columns
        assert "geometry" in gdf.columns
        assert gdf.crs == "EPSG:3857"

    def test_read_admin_parquet(self, test_files):
        """Test reading admin unit Parquet file."""
        gdf = read_geospatial_file(str(test_files["admin_parquet"]))

        assert isinstance(gdf, gpd.GeoDataFrame)
        assert len(gdf) == 6
        assert "ID_admin_unit" in gdf.columns
        assert "geometry" in gdf.columns
        assert gdf.crs == "EPSG:3857"

    def test_read_txt_file_fails(self, test_files):
        """Test that reading TXT file raises error."""
        with pytest.raises(FileNotFoundError, match="Unsupported file type"):
            read_geospatial_file(str(test_files["txt_file"]))

    def test_read_blank_file_fails(self, test_files):
        """Test that reading blank file raises error."""
        with pytest.raises(
            Exception
        ):  # Could be various exceptions depending on file format
            read_geospatial_file(str(test_files["blank_file"]))

    def test_read_nonexistent_file_fails(self):
        """Test that reading nonexistent file raises error."""
        with pytest.raises(Exception):
            read_geospatial_file("nonexistent_file.geojson")

    # ===========================================
    # HAZARD VALIDATION TESTS
    # ===========================================

    def test_validate_hazard_columns_correct_columns(self, hazard_test_data):
        """Test hazard validation passes with correct columns."""
        assert validate_hazard_columns(hazard_test_data) == True

    def test_validate_hazard_columns_missing_id_hazard(self, hazard_test_data):
        """Test hazard validation fails when ID_hazard missing."""
        invalid_data = hazard_test_data.drop(columns=["ID_hazard"])
        assert validate_hazard_columns(invalid_data) == False

    def test_validate_hazard_columns_missing_geometry(self, hazard_test_data):
        """Test hazard validation fails when geometry missing."""
        invalid_data = hazard_test_data.drop(columns=["geometry"])
        assert validate_hazard_columns(invalid_data) == False

    def test_validate_hazard_columns_missing_buffer_dist(self, hazard_test_data):
        """Test hazard validation fails when no buffer_dist columns."""
        invalid_data = hazard_test_data.drop(
            columns=["buffer_dist_500", "buffer_dist_1000"]
        )
        assert validate_hazard_columns(invalid_data) == False

    def test_validate_hazard_columns_one_buffer_dist_ok(self, hazard_test_data):
        """Test hazard validation passes with only one buffer_dist column."""
        single_buffer = hazard_test_data.drop(columns=["buffer_dist_1000"])
        assert validate_hazard_columns(single_buffer) == True

    def test_validate_hazard_columns_empty_dataframe(self):
        """Test hazard validation fails with empty dataframe."""
        empty_gdf = gpd.GeoDataFrame()
        assert validate_hazard_columns(empty_gdf) == False

    # ===========================================
    # ADMIN UNIT VALIDATION TESTS
    # ===========================================

    def test_validate_admin_unit_columns_correct_columns(self, admin_unit_test_data):
        """Test admin unit validation passes with correct columns."""
        assert validate_admin_unit_columns(admin_unit_test_data) == True

    def test_validate_admin_unit_columns_missing_id_admin_unit(
        self, admin_unit_test_data
    ):
        """Test admin unit validation fails when ID_admin_unit missing."""
        invalid_data = admin_unit_test_data.drop(columns=["ID_admin_unit"])
        assert validate_admin_unit_columns(invalid_data) == False

    def test_validate_admin_unit_columns_missing_geometry(self, admin_unit_test_data):
        """Test admin unit validation fails when geometry missing."""
        invalid_data = admin_unit_test_data.drop(columns=["geometry"])
        assert validate_admin_unit_columns(invalid_data) == False

    def test_validate_admin_unit_columns_empty_dataframe(self):
        """Test admin unit validation fails with empty dataframe."""
        empty_gdf = gpd.GeoDataFrame()
        assert validate_admin_unit_columns(empty_gdf) == False

    def test_validate_admin_unit_columns_accepts_perverse_geometries(
        self, admin_unit_test_data
    ):
        """Test admin unit validation passes even with non-polygon geometries."""
        # The validation only checks for required columns, not geometry types
        assert validate_admin_unit_columns(admin_unit_test_data) == True


# Test runner
if __name__ == "__main__":
    pytest.main([__file__])
