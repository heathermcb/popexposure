## `Pop_Exp`: functions to assess the number of people exposed to an environmental hazard

`Pop_Exp` is an open-source Python package designed to help environmental epidemiologists reproducibly assess population-level exposure to environmental hazards based on residential proximity. Functions in the package are designed to be fast, memory-efficient, and easy to use.

### I. Overview

`Pop_Exp` identifies the number of people living near an environmental hazard or set of environmental hazards by overlaying buffered environmental hazard geospatial data with gridded population data, in order to count the number of people that live within the affected area.

`Pop_Exp` can estimate either (a) the number of people living within the buffer distance of each hazard (e.g., the number of people living within 10 km of each individual wildfire disaster burned area in 2018 in California) or (b) the number of people living within the buffer distance of any of the cumulative set of hazards (e.g., the number of people living within 10 km of any wildfire disaster burned area in 2018 in California). These estimates can be broken down by additional spatial units such as census tracts, counties, or ZCTAs. For example, `Pop_Exp` can find the number of people living within 10 km of any wildfire disaster burned area in 2018 by ZCTA, and calculate spatial unit denominators such as the number of residents in each ZCTA. `Pop_Exp` could also be used to assess exposure based on residential proximity to hurricanes, gas wells, emissions sources, or other environmental hazards mapped in geospatial data.

A tutorial on how to use all the functions in `Pop_Exp` is available on GitHub, at [https://github.com/heathermcb/Pop_Exp/demo](https://github.com/heathermcb/Pop_Exp/demo). Please see the tutorial for a detailed explanation on how to use the package functions for exposure assessment.

### II. Available functions

There are three functions available in `Pop_Exp`:

#### 1. `find_num_people_affected`: Estimate the number of people living within a buffer distance of an environmental hazard or set of environmental hazards.

**Inputs:**

- **Environmental hazard geospatial data**: Path to a geospatial data file (GeoJSON or GeoParquet file) with geometries describing environmental hazards (e.g., wildfire disaster boundaries in the US 2015-2020, active oil and gas well coordinates in Texas in 2018, or global data of all tropical cyclone paths in 2020). Hazards can be any kind of geometry except a geometry collection. Columns must include:
  - `ID_climate_hazard`: Unique identifier for each hazard.
  - `geometry`: Geometry of the hazard.
  - `buffer_dist`: Buffer distance to be applied to each hazard. This function buffers each hazard with the corresponding buffer distance passed by the user. If all values in `buffer_dist` are the same, the same buffer distance applied to all hazards. The user can populate the column `buffer_dist` with different values for each hazards, based on attributes of the hazards or hazard size. If all `buffer_dist values` are `0`, no buffer will be applied to the hazards.
- **Population raster data**: Gridded population dataset as a raster file (e.g., GHSL population dataset which provides coverage for the US from 1990-2020 available for download here: [https://human-settlement.emergency.copernicus.eu/ghs_pop.php](https://human-settlement.emergency.copernicus.eu/ghs_pop.php)).
- **Argument 'by_unique_hazard**: This additional argument determines how `find_num_people_affected` will count people living within buffered hazard boundaries.
  - `by_unique_hazard` _(required, True/False)_: When this parameter is set to `True`, `find_num_people_affected` will compute the population affected by each hazard separately. When `True`, if hazards overlap with each other, the same people may be counted as exposed to two or more distinct hazards (e.g., double counted or more). When this parameter is set to `False`, `find_num_people_affected` will compute the total number of people exposed to any hazard in the set of hazards passed to the function. See `Key Features` for more information. There is no default.

**Outputs:**

- Dataframe containing:
  - `ID_climate_hazard`: If `by_unique_hazard` was `True`, this will be a column containing the unique identifiers for each hazard. If `by_unique_hazard` was `False`, this will be a list of unique IDs for hazards that did not overlap with other hazards, and concatenated strings of the IDs of any groups of overlapping hazards.
  - `num_people_affected`: Number of people living within the buffer distance of each hazard or group of hazards, where each row is a single hazard ID or a concatenated list of hazard IDs. Again, if `by_unique_hazard` was `True`, the output data will contain one row for every `ID_climate_hazard`. This means that the rows will be mutually non-exclusive and people may be double counted if they are in the buffered area of two or more different hazards. If `by_unique_hazard` was `False`, the output will contain concatenated `ID_climate_hazard`s wherever hazards or hazard buffers overlapped and people will be counted once if they were in the area of the group of two or more overlapping buffered hazards. See `Key Features` for a detailed explanation.

**Key Features:**

- For overlapping hazard geometries or buffered hazard geometries, the user can choose from two options using the argument `by_unique_hazard`. When this parameter is set to `True`, `find_num_people_affected` estimates the population affected by each hazard separately; if hazards overlap with each other, the same people may be counted as exposed to two or more distinct hazards (i.e., double counted or more). When it is set to `False`, the population affected by overlapping hazards will be combined and people will be counted once if they were in the area of two or more overlapping hazards. The IDs of any overlapping hazards will be concatenated in the output, and the number of people living within the union of those buffered hazards will be returned.
- The user can select the gridded population dataset they want to use based on the population and time period of interest.
- `find_num_people_affected` uses the buffer distances created and passed by the user to buffer the hazard geometries in the best Universal Transverse Mercator projection for each environmental hazard, based on the hazard centroid latitude and longitude, ensuring the most accurate buffered area is created and minimizing distortion from wild map projections.
- `find_num_people_affected` masks the gridded residential population raster with the buffered hazard geometries using partial pixel masking, using the package exactextract. This means if the hazard geometry overlaps with half of a given pixel, when all pixel values within the buffered hazard geometry are summed, half of the given pixel value is added to the sum. This produces the most accurate count of people affected, in contrast to centroid masking or including the entire value of any pixels touched by the hazard geometry.
- Raster masking is done sequentially in exactextract for each individual geometry or set of overlapping geometries in the set of hazards to minimize working memory use and maximize computation speed.

#### 2. `find_num_people_affected_by_geo`: Estimate the number of people living within a buffer distance of an environmental hazard or set of environmental hazards by additional geographies (e.g., census tract, ZCTA).

This function is very similar to `find_num_people_affected`, but returns output by an additional geography (e.g., ZCTAs, counties, census tracts). It provides the number of people living near a buffered hazard or set of buffered hazards in each ZCTA, county, etc. It requires an additional input of a geospatial dataset of additional spatial unit (eg. ZCTA, census tract) geometries.

**Inputs:**

- **Environmental hazard geospatial data**: Path to a geospatial data file (GeoJSON or GeoParquet file) with geometries describing environmental hazards (e.g., wildfire disaster boundaries in the US 2015-2020, active oil and gas well coordinates in Texas in 2018, or global data of all tropical cyclone paths in 2020). Hazards can be any kind of geometry except a geometry collection. Columns must include:
  - `ID_climate_hazard`: Unique identifier for each hazard.
  - `geometry`: Geometry of the hazard.
  - `buffer_dist`: Buffer distance to be applied to each hazard. This function buffers each hazard with the corresponding buffer distance passed by the user. If all values in `buffer_dist` are the same, the same buffer distance applied to all hazards. The user can populate the column `buffer_dist` with different values for each hazards, based on attributes of the hazards or hazard size. If all `buffer_dist values` are `0`, no buffer will be applied to the hazards.
- **Population raster data**: Gridded population dataset as a raster file (e.g., GHSL population dataset which provides coverage for the US from 1990-2020 available for download here: [https://human-settlement.emergency.copernicus.eu/ghs_pop.php](https://human-settlement.emergency.copernicus.eu/ghs_pop.php)).
- **Geographic boundaries**: Path to a geospatial data file containing boundaries an additional spatial unit (e.g., ZCTAs, counties, census tracts). Columns must include:
  - `ID_spatial_unit`: Unique identifier for each geography
  - `geometry`: Geometry of the spatial units
- **Argument 'by_unique_hazard**: This additional argument determines how `find_num_people_affected` will count people living within buffered hazard boundaries.
  - `by_unique_hazard` _(required, True/False)_: When this parameter is set to `True`, `find_num_people_affected` will compute the population affected by each hazard separately. When `True`, if hazards overlap with each other, the same people may be counted as exposed to two or more distinct hazards (e.g., double counted or more). When this parameter is set to `False`, `find_num_people_affected` will compute the total number of people exposed to any hazard in the set of hazards passed to the function. See `Key Features` for more information. There is no default.

**Outputs:**

- Dataframe containing:
  - `ID_climate_hazard`: Unique identifier for each hazard
  - `ID_spatial_unit`: Unique identifier for each additional spatial unit (e.g., ZCTA, county, census tract)
  - `num_people_affected`: Number of people living within the buffer distance of each hazard or group of hazards, where each row is a single hazard ID or a concatenated list of hazard IDs. Again, if `by_unique_hazard` was `True`, the output data will contain one row for every `ID_climate_hazard`. This means that the rows will be mutually non-exclusive and people may be double counted if they are in the buffered area of two or more different hazards. If `by_unique_hazard` was `False`, the output will contain concatenated `ID_climate_hazard`s wherever hazards or hazard buffers overlapped and people will be counted once if they were in the area of the group of two or more overlapping buffered hazards. See `Key Features` for a detailed explanation.

**Key Features:**

This function returns the count of people for each `ID_climate_hazard` - `ID_spatial_unit` combination. For example, a given row may contain the number of people affected by `ID_climate_hazard` 123 in ZCTA 98107. `by_unique_hazard` works the same way here as above in `find_num_people_affected`, as does hazard buffering and partial pixel masking.

#### 3. `find_number_of_people_residing_by_geo`: Estimate the number of people living within additional spatial units such as ZCTAs, counties, or census tracts.

This function is designed to be used with `find_num_people_affected_by_geo` to produce denominators for spatial units.

**Inputs:**

- **Geographic boundaries**: Path to a geospatial data file containing boundaries an additional spatial unit (e.g., ZCTAs, counties, census tracts). Columns must include:
  - `ID_spatial_unit`: Unique identifier for each geography
  - `geometry`: Geometry of the spatial units
- **Population raster data**: Gridded population dataset as a raster file (e.g., GHSL population dataset which provides coverage for the US from 1990-2020 available for download here: [https://human-settlement.emergency.copernicus.eu/ghs_pop.php](https://human-settlement.emergency.copernicus.eu/ghs_pop.php)).

**Outputs:**

- Dataframe containing:
  - `ID_spatial_unit`: Unique identifier for each additional spatial unit (e.g., ZCTA, county, census tract)
  - `num_people_residing`: The number of people living within the boundaries of each spatial unit according to the population raster passed to the function.

### III. Requirements

1. **Python**
   If you do not already have Python, you can install Python at [https://www.python.org/downloads/](https://www.python.org/downloads/).
   We recommend programming in Python with VS Code.
   We've provided a virtual environment containing the requirements of `Pop_Exp` at [https://github.com/heathermcb/PopExp](https://github.com/heathermcb/PopExp).

2. **Inputs**
   You need:

   - A path to a **hazard geospatial data file** (e.g., oil wells, wildfires, floods, etc.) with specific columns.

     columns: `ID_climate_hazard`, `geometry`, `buffer_dist`
     filetype: GeoParquet or GeoJSON

   - A path to a **gridded population dataset**

     format: raster, must contain a CRS

   - A path to **additional geographies geospatial data file** (optional)

     columns: `ID_spatial_unit`, `geometry`
     filetype: GeoParquet or GeoJSON

### IV. How to run

Please see the tutorial for a detailed explanation of how to use all functions in `PopExp`, at [https://github.com/heathermcb/pop_exp/demo](https://github.com/heathermcb/pop_exp/demo).

You can run the function in Python by calling the function with the appropriate arguments. Note that below `hazard_gdf`, `pop_raster`, and `geo_gdf` are the paths to the respective files. Example calls:

```python
python find_num_people_affected(hazard_gdf, pop_raster, buffer_dist=1000)
python find_num_people_affected_by_geo(hazard_gdf, pop_raster, geo_gdf, buffer_dist=1000)
```


### V. Additional must-reads on how Pop_Exp works

Written in plainer language!

**Temporality**: Hazard data is for a specific time period. Maybe you have fracking-related quakes for 2010, or wildfires for 2019. `PopExp` requires you to pick the gridded population raster that you want to use to calculate how many people live near those hazards yourself, so pick one that corresponds to the correct time period. For example, if you have hazard data from 2009-2021, you might not want to use the same population dataset for all of your environmental hazards. You might want to call the function several times on subsets of your data. Maybe you want to call it for each year between 2009 and 2015 using the GHSL population raster from 2010, and then again for each year between 2016 and 2021 with the population raster for 2020. This is up to you to handle.

**Overlapping hazards**: Depending on your dataset, you might have some overlapping hazards. Maybe you are looking at oil wells and you want to know how many people live within 1 km of oil wells in the US. Because there are often multiple wells next to each other, there may be people who live within 1 km of multiple wells. The parameter `by_unique_hazard` allows you to specify how you want to count people. If `by_unique_hazard=False`, the function counts people ONCE if they are within the buffer distance of any hazard, and returns output with overlapping hazards grouped together. It doesn't tell you if people are within the buffer distance of multiple hazards. If `by_unique_hazard=True`, it tells you how many people are within the buffer of each hazard, and double-counts people who are within the buffer of two or more hazards.

This means if you have multiple years or months of data, even if you're using the same population dataset, you might want to do separate function runs for separate years or months. For example, if you have wildfire data from 2015-2020, and you want to know how many people were affected by fires by ZCTA by year, if you throw all the data in this function at once with `by_unique_hazard=False`, if a fire burned ZCTA 10032 in 2015 and in 2020, all you will know from the function output is the count of people who were within the buffer distance of EITHER fire perimeter in that ZCTA. That may not be what you want. The results will not be broken down by year. So you could instead run the function once for each year 2015-2020, to determine how many people were affected by any fire by year.


### VI. Additional details on how these functions work if you're interested

Here are a list of helpers and main functions in the source code of this package. 
Whether these helpers are all called by each main function and in which order changes based on what's being calculated.

1: `prep_geographies`
This function reads in a climate hazard geospatial data file and spatial unit geospatial data file (counties, zcta, etc.) if applicable, in GeoParquet format or GeoJSON format. This file must contains a string column called ID_climate_hazard, a numeric column called 'buffer_dist', and a geometry column, and nothing else. This function makes geometries valid, reprojects to WGS84 projection, and if the data is hazard data, adds a column indicating the best UTM projection to the data frame. 

- Input: dataframe with 3 columns
  1.`ID_climate_hazard` or `ID_spatial_unit`
  2. `buffer_dist`
  3. `geometry`
- Output: dataframe with 4 columns
  1. ID column: `ID_climate_hazard` or `ID_spatial_unit`
  2. original hazard geometry (`geometry`)
  3. buffer distance (`buffer_dist`)
  4. column with best UTM projection for each hazard (not present for spatial unit dataframes) (`utm_projection`)


2: `add_buffered_geom_col`
This function reprojects each hazard into the proper UTM zone and then buffers the provided geometries. The buffer distance is provided in the buffer_dist column. If the buffer distance is 0, the function does not buffer the geometry. After the buffer distance is added, the function reprojects the geometry back to WGS84.

- Input: dataframe with 3 columns
  1. `ID_climate_hazard`
  2. `buffer_dist`
  3. `geometry`
  4. `utm_projection`
- Output: dataframe with 5 columns
  1. `ID_climate_hazard`
  2. `buffer_dist`
  3. `geometry`
  4. `utm_projection`
  5. buffered hazard geometry: `buffered_hazard` 


3: `combine_overlapping_geometries`
This function combines any overlapping hazard geometries into a single geometry, using the GeoPandas function unary_union. This is necessary when calling `find_num_people_affected` or `find_num_people_affected_by_geo` with by_unique_hazard=True. This allows these functions to avoid double counting people who are within the buffered area of two or more hazards. This function retains an ID column that concatenates the IDs of the overlapping hazards.

    Input: prepared data that has been output from prep_data
    Output: dataframe with 3 columns
        climate hazard ID (ID_climate_hazard)
        original hazard geometry (geometry)

4. 

Step 2: `add_buffer_distance_col`
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

























## `Pop_Exp`: functions to assess the number of people exposed to an environmental hazard

`Pop_Exp` is an open-source Python package designed to help environmental epidemiologists reproducibly assess population-level exposure to environmental hazards based on residential proximity. Functions in the package are designed to be fast, memory-efficient, and easy to use.


### I. Overview

`Pop_Exp` identifies the number of people living near an environmental hazard or set of environmental hazards by overlaying buffered environmental hazard geospatial data with gridded population data, in order to count the number of people that live within the affected area. 

`Pop_Exp` can estimate either (a) the number of people living within the buffer distance of each hazard (e.g., the number of people living within 10 km of each individual wildfire disaster burned area in 2018 in California) or (b) the number of people living within the buffer distance of any of the cumulative set of hazards (e.g., the number of people living within 10 km of any wildfire disaster burned area in 2018 in California). These estimates can be broken down by additional spatial units such as census tracts, counties, or ZCTAs. For example, `Pop_Exp` can find the number of people living within 10 km of any wildfire disaster burned area in 2018 by ZCTA, and calculate spatial unit denominators such as the number of residents in each ZCTA. `Pop_Exp` could also be used to assess exposure based on residential proximity to hurricanes, gas wells, emissions sources, or other environmental hazards mapped in geospatial data. 

A tutorial on how to use all the functions in `Pop_Exp` is available on GitHub, at https://github.com/heathermcb/PopExp/demo. Please see the tutorial for a detailed explanation on how to use the package functions for exposure assessment.

### II. Available functions

There are three functions available in `Pop_Exp`:

### 1. `find_num_people_affected`: Estimate the number of people living within an buffer distance of an environmental hazard or set of environmental hazards.

   **Inputs:**

  - **Environmental hazard geospatial data**: Path to a geospatial data file (GeoJSON or GeoParquet file) with geometries describing environmental hazards (e.g., wildfire disaster boundaries in the US 2015-2020, active oil and gas well coordinates in Texas in 2018, or global data of all tropical cyclone paths in 2020). Hazards can be any kind of geometry except a geometry collection. Columns must include:
    - `ID_climate_hazard`: Unique identifier for each hazard.
    - `geometry`: Geometry of the hazard.
    - `buffer_dist`: Buffer distance to be applied to each hazard. This function buffers each hazard with the corresponding buffer distance passed by the user. If all values in `buffer_dist` are the same, the same buffer distance applied to all hazards. The user can populate the column `buffer_dist` with different values for each hazards, based on attributes of the hazards or hazard size. If all `buffer_dist values` are `0`, no buffer will be applied to the hazards. 
  - **Population raster data**: Gridded population dataset as a raster file (e.g., GHSL population dataset which provides coverage for the US from 1990-2020 available for download here: https://human-settlement.emergency.copernicus.eu/ghs_pop.php).
  - **Argument 'by_unique_hazard**: This additional argument determines how `find_num_people_affected` will count people living within buffered hazard boundaries.
  - `by_unique_hazard`_(required, True/False)_: When this parameter is set to `True`, `find_num_people_affected` will compute the population affected by each hazard separately. When `True`, if hazards overlap with each other, the same people may be counted as exposed to two or more distinct hazards (e.g., double counted or more). When this parameter is set to `False`, `find_num_people_affected` will compute the total number of people exposed to any hazard in the set of hazards passed to the function. See `Key Features` for more information. There is no default.

   **Outputs:**

  - Dataframe containing:
    - `ID_climate_hazard`: If `by_unique_hazard` was `True`, this will be a column containing the unique identifiers for each hazard. If `by_unique_hazard` was `False`, this will be a list of unique IDs for hazards that did not overlap with other hazards, and concatenated strings of the IDs of any groups of overlapping hazards.
    - `num_people_affected`: Number of people living within the buffer distance of each hazard or group of hazards, where each row is a single hazard ID or a concatenated list of hazard IDs. Again, if `by_unique_hazard` was `True`, the output data will contain one row for every `ID_climate_hazard`. This means that the rows will be mutually non-exclusive and people may be double counted if they are in the buffered area of two or more different hazards. If `by_unique_hazard` was `False`, the output will contain concatenated `ID_climate_hazard`s wherever hazards or hazard buffers overlapped and people will be counted once if they were in the area of the group of two or more overlapping buffered hazards. See `Key Features` for a detailed explanation.

   **Key Features:**
   
  - For overlapping hazard geometries or buffered hazard geometries, the user can choose from two options using the argument `by_unique_hazard`. When this parameter is set to `True`, `find_num_people_affected` estimates the population affected by each hazard separately; if hazards overlap with each other, the same people may be counted as exposed to two or more distinct hazards (i.e., double counted or more). When it is set to `False`, the population affected by overlapping hazards will be combined and people will be counted once if they were in the area of two or more overlapping hazards. The IDs of any overlapping hazards will be concatenated in the output, and the number of people living within the union of those buffered hazards will be returned.
  - The user can select the gridded population dataset they want to use based on the population and time period of interest.
  -`find_num_people_affected` uses the buffer distances created and passed by the user to buffer the hazard geometries in the best Universal Transverse Mercator projection for each environmental hazard, based on the hazard centroid latitude and longitude, ensuring the most accurate buffered area is created and minimizing distortion from wild map projections. 
  - `find_num_people_affected` masks the gridded residential population raster with the buffered hazard geometries using partial pixel masking, using the package exactextract. This means if the hazard geometry overlaps with half of a given pixel, when all pixel values within the buffered hazard geometry are summed, half of the given pixel value is added to the sum. This produces the most accurate count of people affected, in contrast to centroid masking or including the entire value of any pixels touched by the hazard geometry.
  - Raster masking is done sequentially in exactextract for each individual geometry or set of overlapping geometries in the set of hazards to minimize working memory use and maximize computation speed. 

### 2. `find_num_people_affected_by_geo`: Estimate the number of people living within a buffer distance of an environmental hazard or set of environmental hazards by additional geographies (e.g., census tract, ZCTA).

   This function is very similar to `find_num_people_affected`, but returns output by an additional geography (e.g., ZCTAs, counties, census tracts). It provides the number of people living near a buffered hazard or set of buffered hazards in each ZCTA, county, etc. It requires an additional input of a geospatial dataset of additional spatial unit (eg. ZCTA, census tract) geometries. 

   **Inputs:**

   - **Environmental hazard geospatial data**: Path to a geospatial data file (GeoJSON or GeoParquet file) with geometries describing environmental hazards (e.g., wildfire disaster boundaries in the US 2015-2020, active oil and gas well coordinates in Texas in 2018, or global data of all tropical cyclone paths in 2020). Hazards can be any kind of geometry except a geometry collection. Columns must include:
     - `ID_climate_hazard`: Unique identifier for each hazard.
     - `geometry`: Geometry of the hazard.
     - `buffer_dist`: Buffer distance to be applied to each hazard. This function buffers each hazard with the corresponding buffer distance passed by the user. If all values in `buffer_dist` are the same, the same buffer distance applied to all hazards. The user can populate the column `buffer_dist` with different values for each hazards, based on attributes of the hazards or hazard size. If all `buffer_dist values` are `0`, no buffer will be applied to the hazards. 
   - **Population raster data**: Gridded population dataset as a raster file (e.g., GHSL population dataset which provides coverage for the US from 1990-2020 available for download here: https://human-settlement.emergency.copernicus.eu/ghs_pop.php).

   - **Geographic boundaries**: Path to a geospatial data file containing boundaries an additional spatial unit (e.g., ZCTAs, counties, census tracts). Columns must include:
     - `ID_spatial_unit`: Unique identifier for each geography
     - `geometry`: Geometry of the spatial units
  - **Argument 'by_unique_hazard**: This additional argument determines how `find_num_people_affected` will count people living within buffered hazard boundaries.
  - `by_unique_hazard`_(required, True/False)_: When this parameter is set to `True`, `find_num_people_affected` will compute the population affected by each hazard separately. When `True`, if hazards overlap with each other, the same people may be counted as exposed to two or more distinct hazards (e.g., double counted or more). When this parameter is set to `False`, `find_num_people_affected` will compute the total number of people exposed to any hazard in the set of hazards passed to the function. See `Key Features` for more information. There is no default.

  **Outputs:**

- Dataframe containing:

  - `ID_climate_hazard`: Unique identifier for each hazard
  - `ID_spatial_unit`: Unique identifier for each additional spatial unit (e.g., ZCTA, county, census tract)
     - `num_people_affected`: Number of people living within the buffer distance of each hazard or group of hazards, where each row is a single hazard ID or a concatenated list of hazard IDs. Again, if `by_unique_hazard` was `True`, the output data will contain one row for every `ID_climate_hazard`. This means that the rows will be mutually non-exclusive and people may be double counted if they are in the buffered area of two or more different hazards. If `by_unique_hazard` was `False`, the output will contain concatenated `ID_climate_hazard`s wherever hazards or hazard buffers overlapped and people will be counted once if they were in the area of the group of two or more overlapping buffered hazards. See `Key Features` for a detailed explanation.
  

  **Key Features:**
  
  This function returns the count of people for each `ID_climate_hazard` - `ID_spatial_unit` combination. For example, a given row may contain the number of people affected by `ID_climate_hazard` 123 in ZCTA 98107. `by_unique_hazard` works the same way here as above in `find_num_people_affected`, as does hazard buffering and partial pixel masking. 

### 3. `find_number_of_people_residing_by_geo`: Estimate the number of people living within additional spatial units such as ZCTAs, counties, or census tracts. 

This function is designed to be used with `find_num_people_affected_by_geo` to produce denominators for spatial units. 

   **Inputs:**
   
-  **Geographic boundaries**: Path to a geospatial data file containing boundaries an additional spatial unit (e.g., ZCTAs, counties, census tracts). Columns must include:
     - `ID_spatial_unit`: Unique identifier for each geography
     - `geometry`: Geometry of the spatial units
- **Population raster data**: Gridded population dataset as a raster file (e.g., GHSL population dataset which provides coverage for the US from 1990-2020 available for download here: https://human-settlement.emergency.copernicus.eu/ghs_pop.php).

  **Outputs:**

- Dataframe containing:

  - `ID_spatial_unit`: Unique identifier for each additional spatial unit (e.g., ZCTA, county, census tract)
  - `num_people_residing`: The number of people living within the boundaries of each spatial unit according to the population raster passed to the function.  


### III. Requirements

1. **Python**
   If you do not already have Python, You can install Python at https://www.python.org/downloads/.
   We recommend programming in Python with VS Code.
   We've provided a virtual environment containing the requirements of `PopExp` at https://github.com/heathermcb/PopExp. 

2. **Inputs**
   You need:

   - A path to a **hazard geospatial data file** (eg. oil wells, wildfires, floods, etc.) with specific columns.

   columns: ID_climate_hazard, geometry, buffer_dist
   filetype: GeoParquet or GeoJSON

   - A path to a **gridded population dataset**

   format: raster, must contain a CRS

   - A path to **additional geographies geospatial data file** (optional)

   columns: `ID_spatial_unit`, `geometry`
   filetype: GeoParquet or GeoJSON 


### IV. How to run

Please see the tutorial for a detailed explanation of how to use all functions in `PopExp`, at https://github.com/heathermcb/pop_exp/demo

You can run the function in Python by calling the function with the appropriate arguments. Note that below `hazard_gdf`, `pop_raster`, and `geo_gdf` are the paths to the respective files. Example calls:

  `python find_num_people_affected(hazard_gdf, pop_raster, buffer_dist=1000)`

  `python find_num_people_affected_by_geo(hazard_gdf, pop_raster, geo_gdf, buffer_dist=1000)`

### V. Additional must-reads on how this works

Written in plainer language!

**Temporality**: Hazard data is for a specific time period. Maybe you have fracking-related quakes for 2010, or wildfires for 2019. `PopExp` requires you to pick the gridded population raster that you want to use to calculate how many people live near those hazards yourself, so pick one that corresponds to the correct time period. For example, if you have hazard data from 2009-2021, you might not want to use the same population dataset for all of your environmental hazards. You might want to call the function several times on subsets of your data. Maybe you want to call it for each year between 2009 and 2015 using the GHSL population raster from 2010, and then again for each year between 2016 and 2021 with the population raster for 2020. This is up to you to handle.

**Overlapping hazards**: Depending on your dataset, you might have some overlapping hazards. Maybe you are looking at oil wells and you want to know how many people live within 1 km of oil wells in the US. Because there are often multiple wells next to each other, there may be people who live within 1 km of multiple wells. The parameter `by_unique_hazard` allows you to specify how you want to count people. If `by_unique_hazard=False`, the function counts people ONCE if they are within the buffer distance of any hazard, and returns output with overlapping hazards grouped together. It doesn't tell you if people are within the buffer distance of multiple hazards. If `by_unique_hazard=True`, it tells you how many people are within the buffer of each hazard, and double-counts people who are within the buffer of two or more hazards.

This means if you have multiple years or months of data, even if you're using the same population dataset, you might want to do separate function runs for separate years or months. For example, if you have wildfire data from 2015-2020, and you want to know how many people were affected by fires by ZCTA by year, if you throw all the data in this function at once with `by_unique_hazard=False`, if a fire burned ZCTA 10032 in 2015 and in 2020, all you will know from the function output is the count of people who were within the buffer distance of EITHER fire perimeter in that ZCTA. That may not be what you want. The results will not be broken down by year. So you could instead run the function once for each year 2015-2020, to determine how many people were affected by any fire by year.



------ 

**Main function 1: `find_num_people_affected`**
This function estimates the number of people affected within a buffered area of a hazard. It does this by by overlaying spatial hazard data and gridded population data, using the above helpers.

- Input: path to hazard data, path to population raster, `by_unique_hazard` argument
- Output: dataframe with 2 columns
  1.  climate hazard ID`ID_climate_hazard`
  2.  number of people affected by the hazard (`num_people_affected`)

**Main function 2: `find_num_people_affected_by_geo`**
This function estimates the number of people affected within a buffered area of a hazard by additional geography (e.g., ZCTA, county).

- Input: path to hazard data, path to population raster, path to additional geographies geospatial data file, `by_unique_hazard` argument
- Output: dataframe with 2 columns
  1.  climate hazard ID`ID_climate_hazard`
  2.  ID of the additional geography (`ID_spatial_unit`)
  3.  number of people affected by the hazard (`num_people_affected`)

**Main function 3: `find_number_of_people_residing_by_geo`**
This function estimates the number of people living in each additional spatial unit of a dataframe containing spatial unit boundaries for administrative geographies such as ZCTAs or counties based on a gridded population dataset. 

- Input: path to additional spatial geographies, path to population raster dataset.
- Output: dataframe with 2 columns
  1.  ID of the additional geography (`ID_spatial_unit`)
  2.  Number of people living in that spatial unit (`num_people_residing`)
