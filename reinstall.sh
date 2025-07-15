#!/bin/bash

pip uninstall -y popexposure
pip install .
python -c "import popexposure"
python -c "from popexposure import PopEstimator"