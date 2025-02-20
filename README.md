## PopEx: functions to assess the number of people exposed to an environmental hazard

The goal of PopEx is to determine the number of people residing near environmental hazards such as wildfires, floods, oil wells, or tropical cyclones.

### I. The details

PopEx identifies the number of people affected by an environmental hazard by overlaying environmental hazard geospatial data with gridded population data. The function then counts the number of people that fall within the area of an environmental hazard (or a set of environmental hazards). For example, PopEx can determine the number of people affected by a specific hurricane given spatial information about the area hit by the hurricane and gridded population data. Similarly, PopEx can determine the number of people affected by wildfires in 2024 given spatial information of all areas where wildfires occurred in 2024. As many studies use US census geographies (e.g., census tracts, ZCTAs, etc.), PopEx can also estimate the number of people in the US affected by wildfire disasters from 2015-2020 by ZCTA. _Note: areas affected can be the area of the hazard or a buffered area of the hazard._

To do this, there are two functions available:

1. Estimate the number of people affected by an environmental hazard: `find_num_people_affected`

   **Inputs:**

   - **Environmental hazard geospatial data file**: Path to a geospatial data file (GEOJSON or geoparquet file) containing boundaries of environmental hazards (e.g., wildfire disaster boundaries 2015-2020). Columns must include:
     - `ID_climate_hazard`: Unique identifier for each hazard
     - `geometry`: Polygon geometry of the hazard
     - `buffer_dist`: Buffer distance to be applied to each hazard. Since this is a column in your input data, each row will get whichever buffer is present in that row. If you want the same buffer distance applied to all hazards, you can fill this column in with the same value for all rows. If you want different buffer distances for different hazards, you can fill this column in with different values for each row (e.g., if you want one buffer for small hazards and one for large hazards, you can programmatically fill this column with those two buffer values. Alternatively you could provide a different buffer for every hazard.). If you do not want to use a buffer, you must include this column and fill it with `0`.
   - **Population raster data**: Gridded population dataset (e.g., GHSL population dataset which provides coverage for the US from 1990-2020 available for download here: https://human-settlement.emergency.copernicus.eu/ghs_pop.php)
   - **Function args**: Specification for how people should be counted.

   **Arguments:**

   - `by_unique_hazard`_(required, True/False)_: This parameter is `True` if you want to estimate the population affected by each hazard separately; if hazards overlap with each other, the same people may be counted as exposed to two distinct hazards (e.g., double counted). This parameter should be `False` if you want to union any overlapping hazards and avoid double counting. There is no default.

   **Outputs:**

   - Dataframe containing:
     - `ID_climate_hazard`: Unique identifier for each hazard
     - `num_people_affected`: Number of people affected by each hazard, where each row is a single hazard or a combined hazard. If `by_unique_hazard` was `True`, the output data will contain one row for every `ID_climate_hazard`, and thus the population affected for each hazard. This means that the rows will be mutually non-exclusive and people may be double counted if they are in the buffered area of two or more different hazards. If `by_unique_hazard` was `False`, the output will contain combined `ID_climate_hazard`s whereever hazards or hazard buffers overlapped and people will be counted once if they were in the area of two distinct hazards. See `Key Features` for a detailed explanation.

   **Key Features:**

   - For overlapping hazard geometries or buffered hazard geometries, the user can choose two options:
     1. `by_unique_hazard=True`: Estimates the population affected by each hazard separately; if hazards overlap with each other, the same people may be counted as exposed to two distinct hazards (e.g., double counted)
     2. `by_unique_hazard=False`: Estimates the population affected by overlapping hazards will be combined and people will be counted once if they were in the area of two overlapping hazards. The `ID_climate_hazard` will be a concatenated ID to reflect this.
   - Flexible selection of GHSL population dataset based on the time period of interest
   - Handles multiple hazard types and varying spatial extents

2. Estimate the number of people affected by an environmental hazard by geography (e.g., census tract, ZCTA): `find_num_people_affected_by_geo`

   This function is very similar to `find_num_people_affected`. It differs in that it provides the output by an additional geography (e.g., ZCTAs, counties, census tracts). As such, it has an additional input -- a geospatial data file for the desired geographic aggregation unit.

   **Inputs:**

   - **Environmental hazard geospatial data file**: Path to a geospatial data file containing boundaries of environmental hazards (e.g., wildfire disaster boundaries 2015-2020). Columns must include:
     - `ID_climate_hazard`: Unique identifier for each hazard
     - `geometry`: Polygon geometry of the hazard
     - `buffer_dist`: Buffer distance to be applied to each hazard. Since this is a column in your input data, each row will get whichever buffer is present in that row. If you want the same buffer distance applied to all hazards, you can fill this column in with the same value for all rows. If you want different buffer distances for different hazards, you can fill this column in with different values for each row (e.g., if you want one buffer for small hazards and one for large hazards, you can programmatically fill this column with those two buffer values. Alternatively you could provide a different buffer for every hazard.). If you do not want to use a buffer, you must include this column and fill it with `0`.
   - **Population raster data**: Gridded population dataset (e.g., GHSL population dataset which provides coverage for the US from 1990-2020 available for download here: https://human-settlement.emergency.copernicus.eu/ghs_pop.php)
   - **Geographic boundaries**: Path to a geospatial data file containing boundaries of the desired geographic aggregation unit (e.g., ZCTAs, counties, census tracts). Columns must include:
     - `ID_spatial_unit`: Unique identifier for each geography
     - `geometry`: Polygon geometry of the geography
   - **Function args**: Specification for how people should be counted.

   **Arguments:**

- `by_unique_hazard`_(required, True/False)_: This parameter is `True` if you want to estimate the population affected by each hazard separately; if hazards overlap with each other, the same people may be counted as exposed to two distinct hazards (e.g., double counted). This parameter should be `False` if you want to union any overlapping hazards and avoid double counting. There is no default.

  **Outputs:**

- Dataframe containing:

  - `ID_climate_hazard`: Unique identifier for each hazard
  - `ID_spatial_unit`: Unique identifier for each geography (e.g., ZCTA, county, census tract)
  - `num_people_affected`: Number of people affected by each hazard, where each row is a single hazard or a combined hazard. If `by_unique_hazard` was `True`, the output data will contain one row for every `ID_climate_hazard`, and thus the population affected for each hazard. This means that the rows will be mutually non-exclusive and people may be double counted if they are in the buffered area of two different hazards. If `by_unique_hazard` was `False`, the output will contain combined `ID_climate_hazard`s whereever hazards overlapped and people will be counted once if they were in the area of two distinct hazards. See `Key Features` for a detailed explanation.

- Dataframe containing `ID_climate_hazard`s, ID_spatial_unit (e.g., ZCTA, county, census tract), and the number of people affected by each hazard for that geography. If `by_unique_hazard` was `True`, you will have mutually non-exclusive hazard counts and people may be double counted if they are in the buffered area of two different hazards. If `by_unique_hazard` was `False`, the output data will contain combined `ID_climate_hazard`s whereever hazards overlapped and people will be counted once if they were in the area of two distinct hazards. See `Key Features` for a detailed explanation.

**Key Features:**
This function returns the count of people for each `ID_climate_hazard` - `ID_spatial_unit` combination. For example, a given row may contain the number of people affected by `ID_climate_hazard` 123 in ZCTA 98107.

- For overlapping hazard geometries or buffered hazard geometries, the user can choose two options:
  1.  `by_unique_hazard=True`: Estimates the population affected by each hazard separately; if hazards overlap with each other, the same people may be counted as exposed to two distinct hazards (e.g., double counted)
  2.  `by_unique_hazard=False`: Estimates the population affected by overlapping hazards will be combined and people will be counted once if they were in the area of two overlapping hazards. The `ID_climate_hazard` will be a concatenated ID to reflect this.
- Flexible selection of GHSL population dataset based on the time period of interest
- Handles multiple hazard types and varying spatial extents

### II. Requirements

1. **Python**
   If you do not already have Python, You can install Python at https://www.python.org/downloads/.
   We recommend programming in Python with VS Code.

2. **Inputs**
   You need:

   - **Hazard geospatial data file** (eg. oil wells, wildfires, floods, etc.)

   columns: hazard_id, geometry, buffer_dist
   filetype: parquet or geojson

   - **Gridded population dataset**

   format: raster, must contain a CRS

   - **Additional geographies geospatial data file** (optional)

   columns: hazard_id, geometry
   filetype: parquet or geojson

   note: if you are using census data, you will probably need to clean it before passing it to this function (for example, if I wanted to use census tracts, I'd need to open the raw NHGIS census tract data and rename the GEOID column, and then select that and the geometry column and save that as its own file to pass to this function).

   Note: You need to save both your hazards and additional geometries as parquet or geojson files. You likely have read them in R with the sf package, and they are sf dataframes. If you are using R, you must use a geojson format -- the parquet format will cause issues as it is still under development in R!! If you are using Python, you may use parquet or geojson.

   - In R: You can save them as geojson files by saving those `sf` dataframes using `st_write` in the `sf` package.

### II. How to run

- You can run the function in Python by calling the function with the appropriate arguments. Note that below `hazard_gdf`, `pop_raster`, and `geo_gdf` are the paths to the respective files.

  `python find_num_people_affected(hazard_gdf, pop_raster, buffer_dist=1000)`

  `python find_num_people_affected_by_geo(hazard_gdf, pop_raster, geo_gdf, buffer_dist=1000)`

- Alternatively, you can use a Jupyter notebook to run this.

- Please see the example for how to do this and refer back to the Python tutorial if you need help with this.

**III. Additional must-reads on how this works**

Hazard data is for a specific time period. Maybe you have fracking-related quakes for 2010, or wildfires for 2019. For this function, we're requiring you to pick the gridded population raster that you want to use to calculate how many people live near those hazards yourself, so pick one that corresponds to the correct time period. For example, if you have hazard data from 2009-2021, you might not want to use the same pop dataset for all of your climate hazards. You might want to call the function several times on subsets of your data. Maybe you want to call it for each year between 2009 and 2015 using the GHS pop raster from 2010, and then again for each year between 2016 and 2021 with the pop raster for 2020. This is up to you to handle.

Overlapping hazards: depending on your dataset, you might have some overlapping hazards. Maybe you are looking at oil wells and you want to know how many people live within 1 km of oil wells in the US. Because there are often multiple wells next to each other, there may be people who live within 1 km of multiple wells. The parameter 'by_unique_hazard' allows you to specify how you want to count people. If by_unique_hazard=False, the function counts people ONCE if they are within the buffer distance of any hazard, and returns output with overlapping hazards grouped together. It doesn't tell you if people are within the buffer distance of multiple hazards. If by_unique_hazard=True, it tells you how many people are within the buffer of each hazard, and double-counts people who are within the buffer of two or more hazards.

This means if you have multiple years or months of data, even if you're using the same population dataset, you might want to do separate function runs for separate years or months. For example, if I have wildfire data from 2015-2020, and I want to know how many people were affected by fires by ZCTA by year, if I throw all the data in this function at once with by_unique_hazard=False, if a fire burned ZCTA 10032 in 2015 and in 2020, all I will know is the count of people who were within the buffer distance of EITHER fire perimeter in that ZCTA. That's not what I want. The results will not be broken down by year. So I would run the function once for each year 2015-2020, to know how many people were affected by any fire by year.

**IV. The specifics: more details on how it works if you are interested**

**THIS SECTION IS STILL UNDER DEVELOPMENT**

Step 1: `validate`

TODO: write a validations function to run first and then nothing else happens if it fails

VALIDATIONS: fails if there is more than 2 cols, there is no geometry col, the ID column isnt a string, the geometry isnt a geometry, the polygons are all smaller than the state they are in

Step 2: `prep_geographies`
This function reads in a climate hazard geospatial data file or spatial unit geospatial data file (counties, zcta, etc.) in parquet format that contains a string column with the geom ID and a geography column, but nothing else. This function renames the ID columns consistently, makes geoms valid, and reprojects to albers.

- Input: dataframe should have ONLY the following columns
  `ID_climate_hazard` (can be called whatever you want)
  geometry (must be named `geometry`)
- Output: dataframe is the same but has standardized column names, valid geometries, and is in the albers projection

Step 3: `add_buffer_distance_col`
This function adds a buffer distance to the provided polygons. There are two options for buffer distance, assigned based on whether the hazard area is larger than an area threshold, in square meters.

- Input: dataframe that comes out of `prep_geographies`
- Output: dataframe with 3 columns
  1.  `ID_climate_hazard`
  2.  original hazard geometry (`geometry`)
  3.  buffer distance

Step 4: `add_buffered_geom_col`
This function buffers the original `geometry` using the buffer distance found in step 3.

- Input: dataframe that comes out of `add_buffer_distance_col`
- Output: dataframe with 4 columns
  1.  `ID_climate_hazard`
  2.  original hazard geometry (`geometry`)
  3.  buffer distance
  4.  buffered hazard geometry (`buffered_hazard_geometry`)

Step 5: `add_bounding_box_col`
This function creates a bounding box around the buffered geometry. The purpose is to make it easy to load the population raster only in the bounding box of the hazard for ~_swiftness_~.

- Input: dataframe that comes out of `add_buffered_geom_col`
- Output: dataframe with 5 columns
  1.  `ID_climate_hazard`
  2.  original hazard geometry (`geometry`)
  3.  buffer distance
  4.  buffered hazard geometry (`buffered_hazard_geometry`)
  5.  bounding box around buffered hazard geometry (`bounding_box`)

Step 6: `find_num_people_affected`
This function estimates the number of people affected within a buffered area of a hazard. It does this by **XXX**.

- Input: dataframe that comes out of `add_bounding_box_col`
- Output: dataframe with 6 columns
  1.  `ID_climate_hazard`
  2.  original hazard geometry (`geometry`)
  3.  buffer distance
  4.  buffered hazard geometry (`buffered_hazard_geometry`)
  5.  bounding box around buffered hazard geometry (`bounding_box`)
  6.  number of people affected by the hazard (`num_people_affected`)

Step 7: `find_num_people_affected_by_geo`

Step 8: `find_pop_density`