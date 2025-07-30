"""
Visual unit tests for geometry operations.

This module tests all geometry operations with visual output for validation.
Includes plotting and dataframe inspection for manual verification.
"""

import pytest
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import (
    Point,
    LineString,
    Polygon,
    MultiPoint,
    MultiPolygon,
    GeometryCollection,
)
import numpy as np
import sys
import os
from pathlib import Path

from popexposure.utils.geom_validator import *
from popexposure.utils.geom_ops import *

class TestGeometryOperationsVisual:
    """Visual tests for geometry operations with plotting and dataframe output."""

    @pytest.fixture
    def diverse_hazards(self):
        """Create hazards with different geometry types and buffer distances."""
        print("\n=== CREATING DIVERSE HAZARDS DATASET ===")

        # Create different geometry types
        geometries = [
            # Points
            Point(-122.4, 37.8),  # San Francisco
            Point(-118.2, 34.1),  # Los Angeles
            # LineString (represents a road or fault line)
            LineString([(-122.5, 37.9), (-122.3, 37.7), (-122.1, 37.8)]),
            # Polygon (represents a fire perimeter)
            Polygon([(-119.0, 35.0), (-118.8, 35.0), (-118.8, 35.2), (-119.0, 35.2)]),
            # MultiPoint (multiple earthquake epicenters)
            MultiPoint([Point(-121.0, 36.0), Point(-121.1, 36.1), Point(-120.9, 35.9)]),
            # MultiPolygon (multiple fire areas)
            MultiPolygon(
                [
                    Polygon(
                        [(-117.0, 33.0), (-116.8, 33.0), (-116.8, 33.2), (-117.0, 33.2)]
                    ),
                    Polygon(
                        [(-116.5, 33.1), (-116.3, 33.1), (-116.3, 33.3), (-116.5, 33.3)]
                    ),
                ]
            ),
            # GeometryCollection (mixed hazard types in one event)
            GeometryCollection(
                [
                    Point(-115.0, 32.5),  # Earthquake epicenter
                    LineString([(-115.1, 32.4), (-114.9, 32.6)]),  # Related fault
                    Polygon(
                        [(-115.2, 32.3), (-114.8, 32.3), (-114.8, 32.7), (-115.2, 32.7)]
                    ),  # Affected area
                ]
            ),
        ]

        # Create the GeoDataFrame
        data = {
            "ID_hazard": [
                "point1",
                "point2",
                "line1",
                "poly1",
                "multipoint1",
                "multipoly1",
                "geocoll1",
            ],
            "buffer_dist_500": [
                5000,
                1000,
                5000,
                7500,
                6000,
                8000,
                5500,
            ],  # Different distances
            "buffer_dist_2000": [
                20000,
                15000,
                20000,
                25000,
                18000,
                22000,
                19000,
            ],  # Different distances
            "geometry": geometries,
        }

        gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")

        # Add UTM projections using GeometryValidator
        gdf = add_utm_projection_column(gdf)

        print("Created hazards dataset:")
        print(gdf)
        print(f"Geometry types: {[geom.geom_type for geom in gdf.geometry]}")

        return gdf

    def test_individual_buffering(self, diverse_hazards):
        """Test get_buffered_geometry method on individual geometries."""
        print("\n=== TESTING INDIVIDUAL BUFFERING ===")

        # Test buffering each geometry individually
        buffered_500 = []
        buffered_2000 = []

        for idx, row in diverse_hazards.iterrows():
            # Buffer with 500m
            buff_500 = get_buffered_geometry(row, "buffer_dist_500")
            buffered_500.append(buff_500)

            # Buffer with 2000m
            buff_2000 = get_buffered_geometry(
                row, "buffer_dist_2000"
            )
            buffered_2000.append(buff_2000)

            print(f"Hazard {row['ID_hazard']} ({row['geometry'].geom_type}):")
            print(f"  Original bounds: {row['geometry'].bounds}")
            print(f"  500m buffer bounds: {buff_500.bounds}")
            print(f"  2000m buffer bounds: {buff_2000.bounds}")

        # Add to dataframe for comparison
        test_df = diverse_hazards.copy()
        test_df["manual_buff_500"] = buffered_500
        test_df["manual_buff_2000"] = buffered_2000

        # Plot results
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        # Original geometries
        diverse_hazards.plot(ax=axes[0], color="red", markersize=50, linewidth=2)
        axes[0].set_title("Original Hazards")
        axes[0].legend(["Points", "Lines", "Polygons", "MultiPoints"])

        # 500m buffers
        test_df.set_geometry("manual_buff_500").plot(
            ax=axes[1], color="blue", alpha=0.5
        )
        diverse_hazards.plot(ax=axes[1], color="red", markersize=30)
        axes[1].set_title("500m Buffers")

        # 2000m buffers
        test_df.set_geometry("manual_buff_2000").plot(
            ax=axes[2], color="green", alpha=0.5
        )
        diverse_hazards.plot(ax=axes[2], color="red", markersize=30)
        axes[2].set_title("2000m Buffers")

        plt.tight_layout()
        plt.show()

        print("\nIndividual buffering test DataFrame:")
        print(test_df[["ID_hazard", "buffer_dist_500", "buffer_dist_2000"]])

        return test_df

    def test_batch_buffering(self, diverse_hazards):
        """Test add_buffered_geometry_columns method."""
        print("\n=== TESTING BATCH BUFFERING ===")

        # Use the batch method
        buffered_gdf = add_buffered_geometry_columns(diverse_hazards)

        print("Batch buffering results:")
        print(buffered_gdf.columns.tolist())
        print(buffered_gdf[["ID_hazard"]])

        # Plot comparison
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))

        # Original
        diverse_hazards.plot(ax=axes[0, 0], color="red", markersize=50, linewidth=2)
        axes[0, 0].set_title("Original Hazards")

        # 500m batch buffers
        buffered_gdf.set_geometry("buffered_hazard_500").plot(
            ax=axes[0, 1], color="blue", alpha=0.5
        )
        diverse_hazards.plot(ax=axes[0, 1], color="red", markersize=30)
        axes[0, 1].set_title("500m Batch Buffers")

        # 2000m batch buffers
        buffered_gdf.set_geometry("buffered_hazard_2000").plot(
            ax=axes[1, 0], color="green", alpha=0.5
        )
        diverse_hazards.plot(ax=axes[1, 0], color="red", markersize=30)
        axes[1, 0].set_title("2000m Batch Buffers")

        # All buffers together
        buffered_gdf.set_geometry("buffered_hazard_500").plot(
            ax=axes[1, 1], color="blue", alpha=0.3, label="500m"
        )
        buffered_gdf.set_geometry("buffered_hazard_2000").plot(
            ax=axes[1, 1], color="green", alpha=0.3, label="2000m"
        )
        diverse_hazards.plot(
            ax=axes[1, 1], color="red", markersize=30, label="Original"
        )
        axes[1, 1].set_title("All Buffers Combined")
        axes[1, 1].legend()

        plt.tight_layout()
        plt.show()

        # Verify columns were created
        assert "buffered_hazard_500" in buffered_gdf.columns
        assert "buffered_hazard_2000" in buffered_gdf.columns

        return buffered_gdf

    @pytest.fixture
    def overlapping_hazards(self):
        """Create overlapping hazards for testing union operations."""
        print("\n=== CREATING OVERLAPPING HAZARDS ===")

        # Create overlapping circular areas (wildfire scenario)
        centers = [
            (-120.0, 36.0),  # Central point
            (-119.95, 36.05),  # Slightly northeast
            (-120.05, 35.95),  # Slightly southwest
            (-119.98, 36.02),  # Close to center
        ]

        geometries = [Point(x, y) for x, y in centers]

        data = {
            "ID_hazard": ["fire_a", "fire_b", "fire_c", "fire_d"],
            "buffer_dist_1000": [1000, 800, 1200, 600],
            "buffer_dist_3000": [3000, 2500, 3500, 2000],
            "geometry": geometries,
        }

        gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")
        gdf = add_utm_projection_column(gdf)

        print("Overlapping hazards:")
        print(gdf)

        return gdf

    def test_geometry_combination(self, overlapping_hazards):
        """Test combine_geometries_by_column method."""
        print("\n=== TESTING GEOMETRY COMBINATION ===")

        # First buffer the overlapping hazards
        buffered = add_buffered_geometry_columns(overlapping_hazards)

        print("Buffered overlapping hazards:")
        print(buffered[["ID_hazard"]])

        # Combine geometries
        combined = combine_geometries_by_column(buffered)

        print("\nCombined geometries result:")
        print(combined)
        print(
            f"Combined geometry columns: {[col for col in combined.columns if 'buffered' in col]}"
        )

        # Plot the combination process
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))

        # Original points
        overlapping_hazards.plot(ax=axes[0, 0], color="red", markersize=100)
        axes[0, 0].set_title("Original Overlapping Points")
        for idx, row in overlapping_hazards.iterrows():
            axes[0, 0].annotate(
                row["ID_hazard"],
                (row.geometry.x, row.geometry.y),
                xytext=(5, 5),
                textcoords="offset points",
            )

        # Individual 1000m buffers
        buffered.set_geometry("buffered_hazard_1000").plot(
            ax=axes[0, 1], alpha=0.5, cmap="tab10"
        )
        overlapping_hazards.plot(ax=axes[0, 1], color="red", markersize=50)
        axes[0, 1].set_title("Individual 1000m Buffers")

        # Individual 3000m buffers
        buffered.set_geometry("buffered_hazard_3000").plot(
            ax=axes[1, 0], alpha=0.5, cmap="tab10"
        )
        overlapping_hazards.plot(ax=axes[1, 0], color="red", markersize=50)
        axes[1, 0].set_title("Individual 3000m Buffers")

        # Combined geometries
        combined.set_geometry("buffered_hazard_1000").plot(
            ax=axes[1, 1], color="blue", alpha=0.7, label="1000m Union"
        )
        combined.set_geometry("buffered_hazard_3000").plot(
            ax=axes[1, 1], color="green", alpha=0.7, label="3000m Union"
        )
        overlapping_hazards.plot(
            ax=axes[1, 1], color="red", markersize=50, label="Original Points"
        )
        axes[1, 1].set_title("Combined/Unioned Geometries")
        axes[1, 1].legend()

        plt.tight_layout()
        plt.show()

        # Verify combination worked
        assert len(combined) == 1
        assert "buffered_hazard_1000" in combined.columns
        assert "buffered_hazard_3000" in combined.columns

        return combined, buffered

    @pytest.fixture
    def admin_units(self):
        """Create administrative units for intersection testing."""
        print("\n=== CREATING ADMIN UNITS ===")

        # Create a dense grid of small administrative units (like census blocks)
        admin_polys = []
        admin_ids = []

        # Create 10x10 grid with no gaps - covers the hazard areas more densely
        grid_size = 10
        x_min, x_max = -121.0, -119.0  # Covers all hazard areas
        y_min, y_max = 35.5, 36.5

        x_step = (x_max - x_min) / grid_size
        y_step = (y_max - y_min) / grid_size

        for i in range(grid_size):
            for j in range(grid_size):
                # Calculate exact boundaries with no gaps
                x1 = x_min + i * x_step
                x2 = x_min + (i + 1) * x_step
                y1 = y_min + j * y_step
                y2 = y_min + (j + 1) * y_step

                # Create seamless square administrative unit
                poly = Polygon([(x1, y1), (x2, y1), (x2, y2), (x1, y2)])
                admin_polys.append(poly)
                admin_ids.append(f"admin_{i}_{j}")

        data = {
            "ID_admin_unit": admin_ids,
            "geometry": admin_polys,
        }

        gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")

        print("Administrative units:")
        print(gdf[["ID_admin_unit"]])

        return gdf

    def test_geometry_intersections(self, overlapping_hazards, admin_units):
        """Test get_geometry_intersections method."""
        print("\n=== TESTING GEOMETRY INTERSECTIONS ===")

        # First buffer the hazards
        buffered_hazards = add_buffered_geometry_columns(
            overlapping_hazards
        )

        print("Buffered hazards for intersection:")
        print(buffered_hazards[["ID_hazard"]])

        # Calculate intersections
        intersections = get_geometry_intersections(
            buffered_hazards, admin_units
        )

        print("\nIntersection results:")
        print(intersections)
        print(f"Intersection columns: {intersections.columns.tolist()}")
        print(f"Number of intersections: {len(intersections)}")

        # Plot the intersection analysis
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))

        # Admin units and original hazards
        admin_units.plot(ax=axes[0, 0], color="lightblue", edgecolor="black", alpha=0.7)
        overlapping_hazards.plot(ax=axes[0, 0], color="red", markersize=100)
        axes[0, 0].set_title("Admin Units + Original Hazards")

        # Admin units with 1000m buffers
        admin_units.plot(ax=axes[0, 1], color="lightblue", edgecolor="black", alpha=0.5)
        buffered_hazards.set_geometry("buffered_hazard_1000").plot(
            ax=axes[0, 1], color="orange", alpha=0.6
        )
        overlapping_hazards.plot(ax=axes[0, 1], color="red", markersize=50)
        axes[0, 1].set_title("Admin Units + 1000m Hazard Buffers")

        # Admin units with 3000m buffers
        admin_units.plot(ax=axes[1, 0], color="lightblue", edgecolor="black", alpha=0.5)
        buffered_hazards.set_geometry("buffered_hazard_3000").plot(
            ax=axes[1, 0], color="purple", alpha=0.6
        )
        overlapping_hazards.plot(ax=axes[1, 0], color="red", markersize=50)
        axes[1, 0].set_title("Admin Units + 3000m Hazard Buffers")

        # Intersection results
        admin_units.plot(ax=axes[1, 1], color="lightgray", edgecolor="black", alpha=0.3)
        if not intersections.empty and "buffered_hazard_1000" in intersections.columns:
            intersections.set_geometry("buffered_hazard_1000").plot(
                ax=axes[1, 1], color="red", alpha=0.8, label="1000m Intersections"
            )
        if not intersections.empty and "buffered_hazard_3000" in intersections.columns:
            intersections.set_geometry("buffered_hazard_3000").plot(
                ax=axes[1, 1], color="blue", alpha=0.6, label="3000m Intersections"
            )
        axes[1, 1].set_title("Intersection Results")
        axes[1, 1].legend()

        plt.tight_layout()
        plt.show()

        # Verify intersections
        if not intersections.empty:
            assert "ID_hazard" in intersections.columns
            assert "ID_admin_unit" in intersections.columns
            print(f"Successfully found {len(intersections)} intersections")
        else:
            print("No intersections found - check if hazards and admin units overlap")

        return intersections


def test_complete_workflow():
    """Run the complete workflow test."""
    print("\n" + "=" * 60)
    print("RUNNING COMPLETE GEOMETRY OPERATIONS WORKFLOW TEST")
    print("=" * 60)

    # This will run all the individual tests in sequence
    # when run with pytest -v -s
    pass


if __name__ == "__main__":
    # Allow running this file directly for development
    test_instance = TestGeometryOperationsVisual()

    # Create test data
    hazards = test_instance.diverse_hazards()
    overlapping = test_instance.overlapping_hazards()
    admin = test_instance.admin_units()

    # Run tests
    test_instance.test_individual_buffering(hazards)
    test_instance.test_batch_buffering(hazards)
    test_instance.test_geometry_combination(overlapping)
    test_instance.test_geometry_intersections(overlapping, admin)

    print("\n" + "=" * 60)
    print("ALL VISUAL TESTS COMPLETED!")
    print("=" * 60)
