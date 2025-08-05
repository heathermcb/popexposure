## Quickstart

```python
import pandas as pd
from popexposure import PopEstimator

# List of years and corresponding hazard file paths
years = [2016, 2017, 2018]
hazard_paths = [
    "hazard_2016.geojson",
    "hazard_2017.geojson",
    "hazard_2018.geojson"
]
pop_path = "my_pop_raster.tif"
admin_units_path = "my_admin_units.geojson"

# Instantiate estimator with population data and admin units
pop_est = PopEstimator(
    pop_data=pop_path,
    admin_data=admin_units_path
)

# Find total num ppl residing <= 10km of each hazard in each year
exposed_list = []

for year, hazard_path in zip(years, hazard_paths):
    # Estimate exposed population
    exposed = pop_est.est_exposed_pop(
        hazard_data=hazard_path,
        hazard_specific=False  # set to True if you want per-hazard results
    )
    exposed['year'] = year
    exposed_list.append(exposed)

exposed_df = pd.concat(exposed_list, axis=0)

# Save output
exposed_df.to_parquet("pop_exposed_to_hazards.parquet")
```
