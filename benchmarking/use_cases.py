# NOTE: if using this file with benchmarking/benchmark.py to evaluate
# CPU time consumed by the code, ensure all of the @profile decorators
# are commented out. If instead you are using this with mprof, 
# Uncomment only the @profile decorators. Finally, if line-by-line memory
# profiling is needed, uncomment both the profile import statement below as well
# as the @profile decorators.

# from memory_profiler import profile
from benchmark import benchmark_target
import popexposure as ex
import pandas as pd
import pathlib
import glob

base_dir = pathlib.Path.cwd()
data_dir = f"{base_dir}/docs/tutorials/demo_data"
wildfire_paths = glob.glob(f"{data_dir}/02_interim_data/*fire*")
ghsl_path = f"{data_dir}/01_raw_data/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R5_C8/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R5_C8.tif"
zcta_path = glob.glob(f"{data_dir}/02_interim_data/*zcta*")[0]
start_year = 2016

# @profile
@benchmark_target
def cumulative():
    est = ex.PopEstimator(pop_data=ghsl_path)

    num_exposed_list = []
    for i in range(3):
        num_exposed = est.est_exposed_pop(hazard_data=wildfire_paths[i], hazard_specific=False)
        num_exposed['year'] = start_year + i
        num_exposed_list.append(num_exposed)

    return pd.concat(num_exposed_list, axis=0)

# @profile
@benchmark_target
def cumulative_zcta():
    est = ex.PopEstimator(pop_data=ghsl_path, admin_data=zcta_path)

    num_exposed_list = []
    for i in range(3):
        num_exposed = est.est_exposed_pop(hazard_data=wildfire_paths[i], hazard_specific=False)
        num_exposed['year'] = start_year + i
        num_exposed_list.append(num_exposed)

    return pd.concat(num_exposed_list, axis=0)

# @profile
@benchmark_target
def haz_specific():
    est = ex.PopEstimator(pop_data=ghsl_path)

    num_exposed_list = []
    for i in range(3):
        num_exposed = est.est_exposed_pop(hazard_data=wildfire_paths[i], hazard_specific=True)
        num_exposed['year'] = start_year + i
        num_exposed_list.append(num_exposed)

    return pd.concat(num_exposed_list, axis=0)

# @profile
@benchmark_target
def haz_specific_zcta():
    est = ex.PopEstimator(pop_data=ghsl_path, admin_data=zcta_path)

    num_exposed_list = []
    for i in range(3):
        num_exposed = est.est_exposed_pop(hazard_data=wildfire_paths[i], hazard_specific=True)
        num_exposed['year'] = start_year + i
        num_exposed_list.append(num_exposed)

    return pd.concat(num_exposed_list, axis=0)


# @profile
@benchmark_target
def total_est():
    est = ex.PopEstimator(pop_data=ghsl_path, admin_data=zcta_path)
    return est.est_total_pop()

if __name__ == "__main__":
    cumulative()
    cumulative_zcta()
    haz_specific()
    haz_specific_zcta()
    total_est()