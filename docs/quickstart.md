## Quickstart

```python
import popexposure as ex
import geopandas as gpd
import pandas as pd

# Set paths
my_pop_raster_path = "path/to/my_pop_raster.tif"
my_admin_units_path = "path/to/my_admin_units.geojson"
my_hazard_data_path = "path/to/my_hazards.geojson"

# Read in hazard data and define buffer(s) of interest (in meters)
hazard_gdf = gpd.read_file(my_hazard_data_path)
hazard_gdf["buffer_dist_0km"] = 0
hazard_gdf["buffer_dist_10km"] = 10_000
hazard_gdf["buffer_dist_20km"] = 20_000

# Instantiate estimator
pop_est = ex.PopEstimator(pop_data=my_pop_raster_path, admin_data=my_admin_units_path)

# Find total num ppl residing <= 0, 10, 20km of each hazard in the hazard_gdf
exposed_df = pop_est.est_exposed_pop(
    hazard_data=hazard_gdf,  # can also pass in a filepath here, assuming it has necessary columns
    hazard_specific=False,   # set to True if you need per-hazard results, False for cumulative exposure
)

# Save output
exposed_df.to_parquet("pop_exposed_to_hazards.parquet")
```
