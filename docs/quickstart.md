## Quickstart

```python
import glob
import pandas as pd
from popexposure.find_exposure import PopEstimator

# Instantiate estimator
pop_est = PopEstimator()

# List of years and corresponding hazard file paths
years = [2016, 2017, 2018]
hazard_paths = [
    "hazard_2016.geojson",
    "hazard_2017.geojson",
    "hazard_2018.geojson"
]
pop_path = "my_pop_raster.tif"
admin_units_path = "my_admin_units.geojson"

# Prepare admin units data
admin_units = pop_est.prep_data(admin_units_path, geo_type="admin_unit")

# Find total num ppl residing <= 10km of each hazard in each year
exposed_list = []

for year, hazard_path in zip(years, hazard_paths):
    # Prepare hazard data for this year
    hazards = pop_est.prep_data(hazard_path, geo_type="hazard")
    # Estimate exposed population
    exposed = pop_est.est_exposed_pop(
        pop_path=pop_path,
        hazard_specific=False,  # set to True if you want per-hazard results
        hazards=hazards,
        admin_units=admin_units # leave as None for estimates not by admin unit
    )
    exposed['year'] = year
    exposed_list.append(exposed)

exposed_df = pd.concat(exposed_list, axis=0)

# Save output
exposed_df.to_parquet("pop_exposed_to_hazards.parquet")
```
