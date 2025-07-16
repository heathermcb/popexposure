#!/bin/bash
echo "⏳ Running demo/demo_code/00_download_data.sh..."

# If not running this from the root directory, prompt user to change to the root directory
if [ ! -d "demo/demo_code" ]; then
  echo "❗ Please run this script from the root directory of the project."
  exit 1
fi

# If curl is not installed, prompt the user to install it
if ! command -v curl &> /dev/null; then
  echo "❗ curl is not installed. Please install it to proceed."
  exit 1
fi

# Create demo/demo_data directory if it does not exist
echo "⏳ Step 1️⃣ /4️⃣ : Creating directories:"
echo "└── demo/demo_data"
echo "  ├── 01_raw_data"
echo "  ├── 02_interim_data"
echo "  └── 03_results"
mkdir -p demo/demo_data
mkdir -p demo/demo_data/01_raw_data
mkdir -p demo/demo_data/02_interim_data
mkdir -p demo/demo_data/03_results
printf "✅ Done!\n\n"

# Download datasets
echo "⏳ Step 2️⃣ /4️⃣ : Downloading US wildfire data..."
curl -L -o demo/demo_data/01_raw_data/wfbz_disasters_conus.geojson https://dataverse.harvard.edu/api/access/datafile/10983328
echo "Saved to demo/demo_data/01_raw_data/wfbz_disasters_conus.geojson"
printf "✅ Done! The corresponding description for this dataset can be found at: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/R73R85 \n\n"

echo "⏳ Step 3️⃣ /4️⃣ : Downloading 2020 US ZCTA shapefile..."
curl -L -o demo/demo_data/01_raw_data/tl_2020_us_zcta520.zip https://www2.census.gov/geo/tiger/TIGER2020/ZCTA520/tl_2020_us_zcta520.zip
unzip -o demo/demo_data/01_raw_data/tl_2020_us_zcta520.zip -d demo/demo_data/01_raw_data/tl_2020_us_zcta520
rm demo/demo_data/01_raw_data/tl_2020_us_zcta520.zip
printf "✅ Done! The corresponding description for this dataset can be found at: https://www.census.gov/programs-surveys/geography/guidance/geo-areas/zctas.html \n\n"

echo "⏳ Step 4️⃣ /4️⃣ : Downloading 2020 California residential population data at 100 m resolution..."
curl -L -o demo/demo_data/01_raw_data/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R5_C8.zip https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2020_GLOBE_R2023A_54009_100/V1-0/tiles/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R5_C8.zip
unzip -o demo/demo_data/01_raw_data/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R5_C8.zip -d demo/demo_data/01_raw_data/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R5_C8
rm demo/demo_data/01_raw_data/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R5_C8.zip
printf "✅ Done! The corresponding description for this dataset can be found at: https://human-settlement.emergency.copernicus.eu/download.php?ds=pop \n\n"

echo "✅ Completed demo/demo_code/00_download_data.sh!"