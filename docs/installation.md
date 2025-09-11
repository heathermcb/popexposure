### Installation

#### Install from PyPI (Recommended)

The easiest way to install `popexposure` is via the latest pre-compiled binaries from PyPI with:

```bash
pip install popexposure
```

You can build `popexposure` from source with:

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

#### Conda environment

We provide a [Conda](https://www.anaconda.com/docs/getting-started/miniconda/install) environment YAML [file](https://github.com/heathermcb/popexposure/blob/main/pop_exp.yml) with all the necessary dependencies for the jupyter notebook tutorials. To install and activate the environment, run the below when in the root directory of popexposure repository:

```bash
# When in root directory of popexposure repository
conda env create -f pop_exp.yml
conda activate pop_exp
```

You can then register the above conda environment as a selectable kernel in Jupyter with:

```bash
python -m ipykernel install --user --name pop_exp --display-name "Python (pop_exp)"
```