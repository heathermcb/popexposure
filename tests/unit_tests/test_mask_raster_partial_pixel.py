"""
Unit tests for the RasterExtractor class.

Tests the mask_raster_partial_pixel method with synthetic raster data
and known geometries to verify population extraction calculations.
"""

import pytest
import geopandas as gpd
import pandas as pd
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
import tempfile
import os
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap

from popexposure.utils.geom_validator import *
from popexposure.utils.geom_ops import *
from popexposure.utils.mask_raster_partial_pixel import *


class TestRasterExtractor:
    """Test cases for RasterExtractor class."""

    @pytest.fixture
    def test_raster_file(self):
        """Create a synthetic population raster for testing."""
        # Create a temporary directory and file
        temp_dir = tempfile.mkdtemp()
        raster_path = os.path.join(temp_dir, "test_population.tif")

        # Define raster parameters
        width, height = 100, 100
        # Bounds: longitude from -122 to -121, latitude from 37 to 38
        west, south, east, north = -122.0, 37.0, -121.0, 38.0
        transform = from_bounds(west, south, east, north, width, height)

        # Create population data - higher values in center, lower at edges
        x = np.linspace(0, width - 1, width)
        y = np.linspace(0, height - 1, height)
        X, Y = np.meshgrid(x, y)

        # Create a population distribution (higher in center)
        center_x, center_y = width // 2, height // 2
        pop_data = (
            np.exp(-((X - center_x) ** 2 + (Y - center_y) ** 2) / (width / 4) ** 2)
            * 1000
        )
        pop_data = pop_data.astype(np.float32)

        # Write the raster
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
            dst.write(pop_data, 1)

        # Plot the raster
        plt.figure(figsize=(8, 6))
        plt.imshow(
            pop_data, extent=[west, east, south, north], cmap="viridis", origin="lower"
        )
        plt.colorbar(label="Population Density")
        plt.title("Synthetic Population Raster for Testing")
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.show()

        yield raster_path, pop_data, transform, (west, south, east, north)

        # Cleanup
        if os.path.exists(raster_path):
            os.unlink(raster_path)
        os.rmdir(temp_dir)

    @pytest.fixture
    def test_geometries(self):
        """Create test geometries with diverse geometry types."""
        # Create different geometry types within our raster bounds (-122 to -121, 37 to 38)
        geometries = [
            # Point at center (should have high population)
            Point(-121.5, 37.5),
            # Polygon covering center area
            Polygon(
                [
                    (-121.6, 37.4),
                    (-121.4, 37.4),
                    (-121.4, 37.6),
                    (-121.6, 37.6),
                    (-121.6, 37.4),
                ]
            ),
            # MultiPolygon (two separate areas)
            MultiPolygon(
                [
                    Polygon(
                        [
                            (-121.8, 37.2),
                            (-121.7, 37.2),
                            (-121.7, 37.3),
                            (-121.8, 37.3),
                            (-121.8, 37.2),
                        ]
                    ),
                    Polygon(
                        [
                            (-121.3, 37.7),
                            (-121.2, 37.7),
                            (-121.2, 37.8),
                            (-121.3, 37.8),
                            (-121.3, 37.7),
                        ]
                    ),
                ]
            ),
            # LineString (represents a road or fault line)
            LineString([(-121.9, 37.1), (-121.5, 37.5), (-121.1, 37.9)]),
            # GeometryCollection (mixed types in one feature)
            GeometryCollection(
                [
                    Point(-121.4, 37.3),  # Earthquake epicenter
                    LineString([(-121.45, 37.25), (-121.35, 37.35)]),
                    Polygon(
                        [
                            (-121.5, 37.2),
                            (-121.3, 37.2),
                            (-121.3, 37.4),
                            (-121.5, 37.4),
                            (-121.5, 37.2),
                        ]
                    ),
                ]
            ),
        ]

        # Create GeoDataFrame with buffer distances
        data = {
            "ID_hazard": ["point1", "poly1", "multipoly1", "line1", "geocoll1"],
            "buffer_dist_100": [100, 100, 100, 100, 100],  # Small buffers (meters)
            "buffer_dist_500": [500, 500, 500, 500, 500],  # Medium buffers (meters)
            "geometry": geometries,
        }

        gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")

        # Add UTM projections using GeometryValidator
        gdf = add_utm_projection_column(gdf)
        gdf = add_buffered_geometry_columns(gdf)

        print(f"\nCreated test geometries:")
        for idx, row in gdf.iterrows():
            print(f"  {row['ID_hazard']}: {row['geometry'].geom_type}")

        # Plot the buffered geometries
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))

        # Plot 100m buffered geometries
        ax1 = axes[0]
        if "buffered_hazard_100" in gdf.columns:
            temp_gdf_100 = gpd.GeoDataFrame(
                gdf, geometry="buffered_hazard_100", crs=gdf.crs
            )
            temp_gdf_100.plot(
                ax=ax1, color="blue", alpha=0.4, edgecolor="darkblue", linewidth=2
            )

        # Plot original geometries on top for reference
        gdf.plot(
            ax=ax1,
            color="red",
            alpha=0.8,
            edgecolor="darkred",
            linewidth=2,
            markersize=50,
        )

        # Add labels for each geometry
        for idx, row in gdf.iterrows():
            centroid = row.geometry.centroid
            ax1.annotate(
                f"{row['ID_hazard']}\n({row['geometry'].geom_type})",
                (centroid.x, centroid.y),
                color="white",
                fontsize=9,
                ha="center",
                va="center",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.8),
            )

        ax1.set_title("100m Buffered Geometries (blue) with Originals (red)")
        ax1.set_xlabel("Longitude")
        ax1.set_ylabel("Latitude")
        ax1.grid(True, alpha=0.3)

        # Plot 500m buffered geometries
        ax2 = axes[1]
        if "buffered_hazard_500" in gdf.columns:
            temp_gdf_500 = gpd.GeoDataFrame(
                gdf, geometry="buffered_hazard_500", crs=gdf.crs
            )
            temp_gdf_500.plot(
                ax=ax2, color="orange", alpha=0.4, edgecolor="darkorange", linewidth=2
            )

        # Plot original geometries on top for reference
        gdf.plot(
            ax=ax2,
            color="red",
            alpha=0.8,
            edgecolor="darkred",
            linewidth=2,
            markersize=50,
        )

        # Add labels for each geometry
        for idx, row in gdf.iterrows():
            centroid = row.geometry.centroid
            ax2.annotate(
                f"{row['ID_hazard']}\n({row['geometry'].geom_type})",
                (centroid.x, centroid.y),
                color="white",
                fontsize=9,
                ha="center",
                va="center",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.8),
            )

        ax2.set_title("500m Buffered Geometries (orange) with Originals (red)")
        ax2.set_xlabel("Longitude")
        ax2.set_ylabel("Latitude")
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

        print(gdf)
        print(gdf.columns)        # Create subset for mask_raster_partial_pixel with only needed columns
        mask_raster_gdf = gdf[
            ["ID_hazard", "buffered_hazard_100", "buffered_hazard_500"]
        ].copy()
        
        # Set buffered_hazard_100 as the active geometry column and set CRS
        mask_raster_gdf = gpd.GeoDataFrame(
            mask_raster_gdf, 
            geometry="buffered_hazard_100", 
            crs=gdf.crs
        )
        
        print(
            f"\nColumns for mask_raster_partial_pixel: {mask_raster_gdf.columns.tolist()}"
        )
        print(f"Active geometry column: {mask_raster_gdf.geometry.name}")
        print(f"CRS: {mask_raster_gdf.crs}")

        return mask_raster_gdf

    def test_raster_extraction_basic_functionality(self, test_raster_file, test_geometries):
        """Test basic functionality of raster extraction with diverse geometry types."""
        raster_path, pop_data, transform, bounds = test_raster_file
        
        print(f"\n{'='*60}")
        print("TESTING BASIC RASTER EXTRACTION FUNCTIONALITY")
        print(f"{'='*60}")
        
        # Display input data info
        print(f"\nInput GeoDataFrame:")
        print(f"  Rows: {len(test_geometries)}")
        print(f"  Columns: {test_geometries.columns.tolist()}")
        print(f"  Active geometry: {test_geometries.geometry.name}")
        print(f"  CRS: {test_geometries.crs}")
        
        print(f"\nRaster info:")
        print(f"  Path: {raster_path}")
        print(f"  Bounds: {bounds}")
        print(f"  Shape: {pop_data.shape}")
        print(f"  Population range: {pop_data.min():.1f} to {pop_data.max():.1f}")
        
        # Test the raster extraction
        print(f"\nRunning raster extraction...")
        result = mask_raster_partial_pixel(test_geometries, raster_path)
        
        print(f"\nExtraction results:")
        print(result)
        
        # Basic validation
        print(f"\nValidation checks:")
        
        # Check structure
        assert len(result) == len(test_geometries), f"Expected {len(test_geometries)} rows, got {len(result)}"
        print(f"✓ Correct number of rows: {len(result)}")
        
        # Check required columns exist
        assert "ID_hazard" in result.columns, "Missing ID_hazard column"
        assert "exposed_100" in result.columns, "Missing exposed_100 column" 
        assert "exposed_500" in result.columns, "Missing exposed_500 column"
        print(f"✓ Required columns present: {result.columns.tolist()}")
        
        # Check all values are non-negative
        assert (result["exposed_100"] >= 0).all(), "Found negative values in exposed_100"
        assert (result["exposed_500"] >= 0).all(), "Found negative values in exposed_500"
        print(f"✓ All exposure values are non-negative")
        
        # Check that larger buffers have >= population than smaller buffers
        buffer_comparison_passed = True
        for idx, row in result.iterrows():
            exposed_100 = row["exposed_100"]
            exposed_500 = row["exposed_500"]
            hazard_id = row["ID_hazard"]
            
            if exposed_500 < exposed_100:
                print(f"  ⚠️  {hazard_id}: 500m buffer ({exposed_500:.1f}) < 100m buffer ({exposed_100:.1f})")
                buffer_comparison_passed = False
            else:
                print(f"  ✓ {hazard_id}: 100m={exposed_100:.1f}, 500m={exposed_500:.1f}")
        
        assert buffer_comparison_passed, "Some 500m buffers have less population than 100m buffers"
        print(f"✓ Larger buffers have >= population than smaller buffers")
        
        # Test specific geometry types
        print(f"\nGeometry type analysis:")
        
        # Point should have some exposure (it's at the center of high population)
        point_result = result[result["ID_hazard"] == "point1"]
        point_100 = point_result["exposed_100"].iloc[0]
        point_500 = point_result["exposed_500"].iloc[0]
        assert point_100 > 0, "Point at center should have positive exposure"
        assert point_500 > point_100, "Point 500m buffer should have more exposure than 100m"
        print(f"  Point (center): 100m={point_100:.1f}, 500m={point_500:.1f} ✓")
        
        # Polygon should have some exposure 
        poly_result = result[result["ID_hazard"] == "poly1"]
        poly_100 = poly_result["exposed_100"].iloc[0]
        poly_500 = poly_result["exposed_500"].iloc[0]
        assert poly_100 > 0, "Polygon should have positive exposure"
        print(f"  Polygon: 100m={poly_100:.1f}, 500m={poly_500:.1f} ✓")
        
        # MultiPolygon (edge areas) might have lower exposure
        multi_result = result[result["ID_hazard"] == "multipoly1"]
        multi_100 = multi_result["exposed_100"].iloc[0]
        multi_500 = multi_result["exposed_500"].iloc[0]
        print(f"  MultiPolygon: 100m={multi_100:.1f}, 500m={multi_500:.1f} ✓")
        
        # LineString should benefit significantly from buffering
        line_result = result[result["ID_hazard"] == "line1"]
        line_100 = line_result["exposed_100"].iloc[0]
        line_500 = line_result["exposed_500"].iloc[0]
        print(f"  LineString: 100m={line_100:.1f}, 500m={line_500:.1f} ✓")
        
        # GeometryCollection should have some exposure
        gc_result = result[result["ID_hazard"] == "geocoll1"]
        gc_100 = gc_result["exposed_100"].iloc[0]
        gc_500 = gc_result["exposed_500"].iloc[0]
        print(f"  GeometryCollection: 100m={gc_100:.1f}, 500m={gc_500:.1f} ✓")
        
        # Summary statistics
        print(f"\nSummary statistics:")
        print(f"  Total 100m exposure: {result['exposed_100'].sum():.1f}")
        print(f"  Total 500m exposure: {result['exposed_500'].sum():.1f}")
        print(f"  Average 100m exposure: {result['exposed_100'].mean():.1f}")
        print(f"  Average 500m exposure: {result['exposed_500'].mean():.1f}")
        print(f"  Max 100m exposure: {result['exposed_100'].max():.1f}")
        print(f"  Max 500m exposure: {result['exposed_500'].max():.1f}")
        
        # Validate that we're extracting reasonable amounts of population
        total_raster_pop = pop_data.sum()
        total_extracted_500 = result['exposed_500'].sum()
        extraction_ratio = total_extracted_500 / total_raster_pop
        
        print(f"\nPopulation extraction validation:")
        print(f"  Total raster population: {total_raster_pop:.1f}")
        print(f"  Total extracted (500m): {total_extracted_500:.1f}")
        print(f"  Extraction ratio: {extraction_ratio:.3f}")
        
        # Should extract some but not all population (geometries don't cover entire raster)
        assert total_extracted_500 > 0, "Should extract some population"
        assert total_extracted_500 <= total_raster_pop * 1.1, "Shouldn't extract more than total raster population (allowing small margin)"
        
        print(f"\n{'='*60}")
        print("✅ ALL BASIC RASTER EXTRACTION TESTS PASSED!")
        print(f"{'='*60}")
        
        return result
