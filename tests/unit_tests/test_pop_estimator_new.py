"""
Unit tests for the new PopEstimator interface.
"""

import pytest
import geopandas as gpd
import numpy as np
import os
import tempfile
import rasterio
from shapely.geometry import Point, Polygon
from popexposure.pop_estimator import PopEstimator

@pytest.fixture
def temp_dir():
    d = tempfile.mkdtemp()
    yield d
    # Clean up
    try:
        import shutil
        shutil.rmtree(d)
    except Exception:
        pass

@pytest.fixture
def standard_raster(temp_dir):
    """Create a 100x100 raster with all 1s covering -122 to -121, 37 to 38."""
    raster_path = os.path.join(temp_dir, "raster_ones.tif")
    width, height = 100, 100
    west, south, east, north = -122.0, 37.0, -121.0, 38.0
    transform = rasterio.transform.from_bounds(west, south, east, north, width, height)
    raster_data = np.ones((height, width), dtype=np.float32)
    with rasterio.open(
        raster_path, "w", driver="GTiff", height=height, width=width, count=1,
        dtype=np.float32, crs="EPSG:4326", transform=transform, nodata=0
    ) as dst:
        dst.write(raster_data, 1)
    return raster_path

@pytest.fixture
def admin_units():
    """Create 2 simple admin polygons."""
    return gpd.GeoDataFrame({
        "ID_admin_unit": ["A", "B"],
        "geometry": [
            Polygon([(-122, 37), (-121.5, 37), (-121.5, 38), (-122, 38), (-122, 37)]),
            Polygon([(-121.5, 37), (-121, 37), (-121, 38), (-121.5, 38), (-121.5, 37)])
        ]
    }, crs="EPSG:4326")

@pytest.fixture
def hazard_data():
    """Create 2 point hazards with 500m and 1000m buffers."""
    return gpd.GeoDataFrame({
        "ID_hazard": ["haz1", "haz2"],
        "buffer_dist_500": [500, 500],
        "buffer_dist_1000": [1000, 1000],
        "geometry": [Point(-121.75, 37.5), Point(-121.25, 37.5)]
    }, crs="EPSG:4326")

@pytest.mark.parametrize("hazard_specific", [True, False])
def test_est_exposed_pop_basic(standard_raster, hazard_data, hazard_specific):
    """Test est_exposed_pop with and without hazard_specific, no admin units."""
    estimator = PopEstimator(pop_data=standard_raster)
    result = estimator.est_exposed_pop(hazard_data, hazard_specific=hazard_specific)
    assert isinstance(result, (gpd.GeoDataFrame, gpd.pd.DataFrame))
    # Should have exposed columns
    assert any(col.startswith("exposed_") for col in result.columns)
    # Should have ID_hazard if hazard_specific
    if hazard_specific:
        assert "ID_hazard" in result.columns
    else:
        assert "ID_hazard" in result.columns or result.shape[0] == 1

@pytest.mark.parametrize("hazard_specific", [True, False])
def test_est_exposed_pop_with_admin(standard_raster, hazard_data, admin_units, hazard_specific):
    """Test est_exposed_pop with admin units."""
    estimator = PopEstimator(pop_data=standard_raster, admin_data=admin_units)
    result = estimator.est_exposed_pop(hazard_data, hazard_specific=hazard_specific)
    assert "ID_admin_unit" in result.columns
    assert any(col.startswith("exposed_") for col in result.columns)
    # Should have one row per admin unit (or per hazard-admin pair)
    if hazard_specific:
        assert set(result["ID_admin_unit"]).issubset({"A", "B"})
    else:
        assert set(result["ID_admin_unit"]).issubset({"A", "B"})

@pytest.mark.parametrize("stat", ["sum", "mean"])
def test_est_total_pop(standard_raster, admin_units, stat):
    """Test est_total_pop returns correct columns and row count."""
    estimator = PopEstimator(pop_data=standard_raster, admin_data=admin_units)
    result = estimator.est_total_pop(stat=stat)
    assert set(result.columns) == {"ID_admin_unit", "population"}
    assert len(result) == 2


def test_error_on_missing_admin_for_total_pop(standard_raster):
    """est_total_pop should raise if admin_data is missing."""
    estimator = PopEstimator(pop_data=standard_raster)
    with pytest.raises(Exception):
        estimator.est_total_pop()


def test_error_on_invalid_hazard_data(standard_raster, admin_units):
    """est_exposed_pop should raise on invalid hazard data."""
    estimator = PopEstimator(pop_data=standard_raster, admin_data=admin_units)
    # Pass a DataFrame missing required columns
    bad_hazards = gpd.GeoDataFrame({"geometry": [Point(-121.5, 37.5)]}, crs="EPSG:4326")
    with pytest.raises(Exception):
        estimator.est_exposed_pop(bad_hazards)
