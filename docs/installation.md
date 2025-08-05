### Installation

#### Install from PyPI (Recommended)

The easiest way to install `popexposure` is via the latest pre-compiled binaries from PyPI with:

```bash
pip install popexposure
```

You can build `popexposure` from source as you would any other Python package with:

```bash
git clone https://github.com/heathermcb/popexposure
cd popexposure
python -m pip install .
```

#### Requirements

`popexposure` requires Python 3.11 or later and depends on several geospatial libraries:

- pandas (any version)
- geopandas == 1.0.1 (pinned version)
- scipy == 1.14.0 (pinned version)
- numpy == 2.0.1 (pinned version)
- matplotlib (any version)
- shapely (any version)
- affine (any version)
- rasterio (any version)
- tqdm (any version)
- pyarrow (any version)
- exactextract (any version)
- pyproj (any version)
