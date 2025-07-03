"""
Unit tests for the GeometryValidator class.

Tests all methods in the GeometryValidator class including geometry cleaning,
CRS operations, and UTM projection calculations.
"""

import pytest
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, Polygon, GeometryCollection
import numpy as np

import sys
from pathlib import Path

# Add the directory containing the module to Python path
module_path = Path("/Volumes/squirrel-utopia/popexposure/src/popexposure")
sys.path.insert(0, str(module_path))

# Import the module directly
import geometry_validator
from geometry_validator import GeometryValidator


class TestGeometryValidator:
    """Test cases for GeometryValidator class."""

    def test_remove_missing_geometries(self):
        """Test removal of rows with null or empty geometries."""
        # Create hazard test data with missing geometries (reusing data structure from data_loader tests)
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
            # Missing/empty geometries (the key test cases)
            None,  # Null geometry
            Point(0, 0).buffer(0).boundary,  # Empty LineString geometry
            Point(0, 0)
            .buffer(0)
            .boundary.buffer(0),  # Empty geometry created by operations
            # Additional perverse case
            Polygon(
                [
                    (-8234000, 4967000),
                    (-8232000, 4969000),
                    (-8232000, 4967000),
                    (-8234000, 4969000),
                    (-8234000, 4967000),
                ]
            ),  # Self-intersecting bowtie (invalid but not empty - should be kept)
        ]

        data = {
            "ID_hazard": ["h1", "h2", "h3", "h4", "h5", "h6"],
            "buffer_dist_500": [500, 750, 0, 0, 0, 1200],
            "buffer_dist_1000": [1000, 1500, 0, 0, 0, 2400],
            "geometry": geometries,
        }
        gdf = gpd.GeoDataFrame(data, crs="EPSG:3857")

        original_length = len(gdf)
        cleaned = GeometryValidator.remove_missing_geometries(gdf)

        # Should have fewer rows after removing missing/empty geometries
        assert len(cleaned) < original_length
        # All remaining geometries should be non-null and non-empty
        assert cleaned["geometry"].notnull().all()
        assert (~cleaned["geometry"].is_empty).all()
        # Should keep the valid geometries (h1, h2, and h6 which is invalid but not missing/empty)
        assert len(cleaned) == 3

    def test_clean_geometries_with_invalid_geoms(self):
        """Test cleaning of invalid geometries makes them all valid."""
        # Create hazard test data with invalid geometries (reusing same structure)
        geometries = [
            # Valid cases
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
            # Invalid geometries to be cleaned
            Polygon(
                [
                    (-8234000, 4967000),
                    (-8232000, 4969000),
                    (-8232000, 4967000),
                    (-8234000, 4969000),
                    (-8234000, 4967000),
                ]
            ),  # Self-intersecting bowtie
            Polygon(
                [(0, 0), (2, 2), (2, 0), (0, 2), (0, 0)]
            ),  # Another self-intersecting
            Polygon(
                [(0, 0), (1, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
            ),  # Duplicate consecutive points
            Polygon([(0, 0), (1, 0), (0, 0)]),  # Collinear points - degenerates to line
            GeometryCollection(
                [
                    Point(100, 100),
                    Polygon(
                        [(100, 101), (101, 101), (101, 102), (100, 102), (100, 101)]
                    ),
                ]
            ),  # GeometryCollection with mixed valid geometries
            GeometryCollection(
                [
                    Point(200, 200),
                    Polygon(
                        [(0, 0), (2, 2), (2, 0), (0, 2), (0, 0)]
                    ),  # Self-intersecting bowtie
                    Polygon(
                        [(0, 0), (1, 0), (0, 0)]
                    ),  # Collinear points - degenerates to line
                ]
            ),  # GeometryCollection with invalid geometries
        ]

        data = {
            "ID_hazard": ["h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8"],
            "buffer_dist_500": [500, 750, 1200, 800, 600, 600, 700, 900],
            "buffer_dist_1000": [1000, 1500, 2400, 1600, 1200, 1200, 1400, 1800],
            "geometry": geometries,
        }
        gdf2 = gpd.GeoDataFrame(data, crs="EPSG:3857")

        assert not gdf2.geometry.is_valid.all()

        # Clean the geometries
        cleaned_gdf = GeometryValidator.clean_geometries(gdf2)

        # All geometries should be valid after cleaning
        assert cleaned_gdf.geometry.is_valid.all()
        # Should have non-empty geometries
        assert (~cleaned_gdf.geometry.is_empty).all()
        # Should keep most geometries since they can be fixed, but expect some losses
        # from collinear polygons and potentially problematic GeometryCollections
        assert len(cleaned_gdf) >= len(gdf2) - 3

    def test_reproject_to_wgs84(self):
        """Test reprojection from non-WGS84 to WGS84."""
        # Create test data in Web Mercator (EPSG:3857)
        web_mercator_point = Point(-8238310, 4969803)  # Roughly NYC in Web Mercator
        data = {"ID": ["test_point"], "geometry": [web_mercator_point]}
        gdf = gpd.GeoDataFrame(data, crs="EPSG:3857")  # Web Mercator

        # Verify starting CRS is not WGS84
        assert gdf.crs != "EPSG:4326"

        # Reproject to WGS84
        wgs84_gdf = GeometryValidator.reproject_to_wgs84(gdf)

        # Verify CRS is now WGS84
        assert str(wgs84_gdf.crs) == "EPSG:4326"

        # Test with data already in WGS84 (should return unchanged)
        wgs84_data = gpd.GeoDataFrame(
            {"ID": ["already_wgs84"], "geometry": [Point(-74, 40.7)]}, crs="EPSG:4326"
        )
        result = GeometryValidator.reproject_to_wgs84(wgs84_data)
        assert str(result.crs) == "EPSG:4326"

    def test_utm_zones_different_locations(self):
        """Test that geometries in different locations get different UTM zones."""
        # Create geometries in different UTM zones (mix of points and polygons)
        locations = [
            Point(-74.0060, 40.7128),  # NYC - UTM Zone 18N
            Point(-122.4194, 37.7749),  # San Francisco - UTM Zone 10N
            Polygon(
                [
                    (2.3522, 48.8566),
                    (2.36, 48.86),
                    (2.36, 48.85),
                    (2.35, 48.85),
                    (2.3522, 48.8566),
                ]
            ),  # Paris polygon - UTM Zone 31N
            Polygon(
                [
                    (151.2093, -33.8688),
                    (151.22, -33.87),
                    (151.22, -33.86),
                    (151.20, -33.86),
                    (151.2093, -33.8688),
                ]
            ),  # Sydney polygon - UTM Zone 56S
            Point(-58.3816, -34.6037),  # Buenos Aires - UTM Zone 21S
            Polygon(
                [
                    (72.8777, 19.0760),
                    (72.9, 19.08),
                    (72.9, 19.07),
                    (72.87, 19.07),
                    (72.8777, 19.0760),
                ]
            ),  # Mumbai polygon - UTM Zone 43N
            Point(139.6917, 35.6895),  # Tokyo - UTM Zone 54N
        ]

        data = {
            "ID_hazard": [
                "nyc",
                "sf",
                "paris",
                "sydney",
                "buenos_aires",
                "london",
                "tokyo",
            ],
            "geometry": locations,
        }
        gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")

        # Add UTM projection columns
        result = GeometryValidator.add_utm_projection_column(gdf)

        # Get unique UTM projections
        utm_projections = result["utm_projection"].unique()

        # Should have different UTM zones for different locations
        assert len(utm_projections) == 7

        # Verify specific expected UTM zones
        nyc_utm = result[result["ID_hazard"] == "nyc"]["utm_projection"].iloc[0]
        sydney_utm = result[result["ID_hazard"] == "sydney"]["utm_projection"].iloc[0]
        paris_utm = result[result["ID_hazard"] == "paris"]["utm_projection"].iloc[0]
        london_utm = result[result["ID_hazard"] == "london"]["utm_projection"].iloc[0]

        # NYC should be in northern hemisphere, Sydney in southern
        assert "326" in nyc_utm  # Northern hemisphere
        assert "327" in sydney_utm  # Southern hemisphere
        assert "326" in paris_utm  # Northern hemisphere (Paris polygon)

    def test_add_utm_projection_column(self):
        """Test that add_utm_projection_column adds the expected columns."""
        # Create test data
        data = {
            "ID_hazard": ["test1", "test2"],
            "geometry": [Point(-74, 40.7), Point(151, -33.9)],
        }
        gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")

        # Add UTM projection columns
        result = GeometryValidator.add_utm_projection_column(gdf)

        # Check that new columns were added
        expected_columns = ["ID_hazard", "geometry", "utm_projection"]
        for col in expected_columns:
            assert col in result.columns

        # Check that utm_projection column contains valid EPSG codes
        utm_projections = result["utm_projection"]
        for utm_code in utm_projections:
            assert utm_code.startswith("EPSG:")
            assert utm_code[5:].isdigit()  # Should be numeric after 'EPSG:'

    def test_get_best_utm_projection_specific_locations(self):
        """Test UTM projection calculation for specific known locations."""
        # Test NYC coordinates
        nyc_utm = GeometryValidator.get_best_utm_projection(40.7128, -74.0060)
        assert nyc_utm == "EPSG:32618"  # UTM Zone 18N

        # Test Sydney coordinates
        sydney_utm = GeometryValidator.get_best_utm_projection(-33.8688, 151.2093)
        assert sydney_utm == "EPSG:32756"  # UTM Zone 56S

        # Test London coordinates
        london_utm = GeometryValidator.get_best_utm_projection(51.5074, -0.1278)
        assert london_utm == "EPSG:32630"  # UTM Zone 30N
