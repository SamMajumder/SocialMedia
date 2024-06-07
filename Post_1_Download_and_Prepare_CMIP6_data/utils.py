# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 07:45:59 2024

@author: Dr. Sambadi Majumder 
"""

### retrieval of the data 
import boto3
from botocore import UNSIGNED
from botocore.client import Config
import s3fs


## preparation of the data 
import osmnx as ox
import rioxarray
import h5netcdf
import xarray as xr
import netCDF4 as nc  


## to work in parallel
from concurrent.futures import ThreadPoolExecutor, as_completed 


## to handle file systems 
import os


########### 
## 


######## 
## Function to extract the bounds based on country or continent 
#############

def extract_bounds_osmnx(place = "India"):
    
    # Get a geodataframe related to the place
    area = ox.geocode_to_gdf(place)  
    
    ### isolate the geometry
    area = area[["geometry"]] 
    
    ## project it to a epsg 4326
    area = area.to_crs("EPSG: 4326") 
    
    ### extracting the bounds
    bounds = area.total_bounds  
    
    ### create a tuple
    bounds = tuple(bounds)

    
    return(bounds)


############ 
### 

#########calculate rolling average for each month
## 

def prepare_cmip6_netcdf_monthly(start_year, end_year, model, scenario, variables, bounds=None, num_workers=None, default_crs=None, target_dataset_crs=None):
    start_year = int(start_year)
    end_year = int(end_year)
    num_workers = int(num_workers) if num_workers is not None else 1
    
    valid_models = ["ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CESM2-WACCM",
                    "CESM2", "CMCC-CM2-SR5", "CMCC-ESM2", "CNRM-CM6-1",
                    "CNRM-ESM2-1", "EC-Earth3-Veg-LR", "EC-Earth3", "FGOALS-g3",
                    "GFDL-CM4", "GFDL-CM4_gr1", "GFDL-ESM4", "GISS-E2-1-G",
                    "HadGEM3-GC31-LL", "HadGEM3-GC31-MM", "IITM-ESM", "INM-CM4-8",
                    "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "KIOST-ESM",
                    "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
                    "MRI-ESM2-0", "NESM3", "NorESM2-LM", "NorESM2-MM", "TaiESM1",
                    "UKESM1-0-LL"]
    
    if model not in valid_models:
        raise ValueError(f"Invalid model '{model}'. Choose from {valid_models}.")

    valid_scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
    if scenario not in valid_scenarios:
        raise ValueError(f"Invalid scenario '{scenario}'. Choose from {valid_scenarios}.")

    fs = s3fs.S3FileSystem(anon=True)

    ensemble_mapping = {
        "CESM2": "r4i1p1f1",
        "CNRM-CM6-1": "r1i1p1f2",
        "GISS-E2-1-G": "r1i1p1f2",
        "MIROC-ES2L": "r1i1p1f2",
        "UKESM1-0-LL": "r1i1p1f2",
        "FGOALS-g3": "r3i1p1f1",
        "HadGEM3-GC31-LL": "r1i1p1f3",
        "HadGEM3-GC31-MM": "r1i1p1f3"
    }
    ensemble = ensemble_mapping.get(model, "r1i1p1f1")

    def get_base_path(variable):
        return f's3://nex-gddp-cmip6/NEX-GDDP-CMIP6/{model}/{scenario}/{ensemble}/{variable}/'

    def process_file(file_path, bounds, variable, fs, default_crs, target_dataset_crs):
        try:
            with fs.open(file_path, mode='rb') as f:
                ds = xr.open_dataset(f, engine='h5netcdf')
                ds.load()
                
                # Check if CRS is set and assign default CRS if not set
                if not ds.rio.crs:
                    if default_crs:
                        ds.rio.write_crs(default_crs, inplace=True)
                
                # Re-project to target CRS if specified and different
                if target_dataset_crs and ds.rio.crs != target_dataset_crs:
                    ds = ds.rio.reproject(target_dataset_crs)
                
                monthly_nc = ds[variable].resample(time='1ME').mean().astype('float32')
                monthly_nc_ds = monthly_nc.to_dataset()
                
                monthly_nc_ds.attrs = ds.attrs
                monthly_nc_ds[variable].attrs = ds[variable].attrs
                
                if bounds:
                    monthly_nc_ds = monthly_nc_ds.rio.clip_box(*bounds)
                
            return monthly_nc_ds

        except Exception as e:
            print(f"Failed to process file {file_path}: {e}")
            return None

    combined_nc_all_vars = None
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        for variable in variables:
            base_path = get_base_path(variable)
            futures = []
            for year in range(start_year, end_year + 1):
                file_path = f"{base_path}{variable}_day_{model}_{scenario}_{ensemble}_gn_{year}.nc"
                if "_v1.1" not in file_path:
                    future = executor.submit(process_file, file_path, bounds, variable, fs, default_crs, target_dataset_crs)
                    futures.append(future)

            processed_datasets = [future.result() for future in as_completed(futures) if future.result() is not None]
            if processed_datasets:
                combined_nc = xr.concat(processed_datasets, dim="time")
                combined_nc['month'] = combined_nc['time'].dt.month
                combined_nc['year'] = combined_nc['time'].dt.year
                combined_nc = combined_nc.set_coords('month')  # Set 'month' as a coordinate
                combined_nc_mean = combined_nc.groupby('month', squeeze=False).mean(dim='time')

                if combined_nc_all_vars is None:
                    combined_nc_all_vars = combined_nc_mean
                else:
                    combined_nc_all_vars = xr.merge([combined_nc_all_vars, combined_nc_mean])
                print(f"All datasets for variable {variable} processed successfully.")
            else:
                print(f"No datasets were processed for variable {variable}.")
        
    return combined_nc_all_vars


############# 

def process_multiple_scenarios_monthly(start_year, end_year, model, variables, scenarios, bounds=None, num_workers=4, log_func=None, default_crs=None, target_dataset_crs=None):
    results = {}
    for scenario in scenarios:
        key = f"{model}_{scenario}"
        message = f"Starting processing for: {key}"
        if log_func:
            log_func(message)
        else:
            print(message)
        try:
            dataset = prepare_cmip6_netcdf_monthly(start_year, end_year, model, scenario, variables, bounds, num_workers, default_crs, target_dataset_crs)
            results[key] = dataset
            if dataset is not None:
                message = f"Data processed for: {key}"
                if log_func:
                    log_func(message)
                else:
                    print(message)
            else:
                message = f"No data returned for: {key}"
                if log_func:
                    log_func(message)
                else:
                    print(message)
        except Exception as e:
            message = f"Exception while processing {key}: {str(e)}"
            if log_func:
                log_func(message)
            else:
                print(message)
            results[key] = None
    return results

###########

def process_multiple_models_monthly(start_year, end_year, models, variables, scenarios, bounds=None, num_workers=4, log_func=None, default_crs=None, target_dataset_crs=None):
    all_results = {}
    for model in models:
        message = f"Starting processing for model: {model}"
        if log_func:
            log_func(message)
        else:
            print(message)
        
        results = process_multiple_scenarios_monthly(start_year, end_year, model, variables, scenarios, bounds, num_workers, log_func, default_crs, target_dataset_crs)
        all_results.update(results)
        
        message = f"Finished processing for model: {model}"
        if log_func:
            log_func(message)
        else:
            print(message)
    
    return all_results



##################### export each dataset by netcdf based on key 

def export_datasets_to_netcdf(dataset_dict, output_folder):
    os.makedirs(output_folder, exist_ok=True)

    for dataset_name, dataset in dataset_dict.items():
        if dataset is not None:
            output_path = os.path.join(output_folder, f"{dataset_name}.nc")
            dataset.to_netcdf(output_path)
            print(f"Saved dataset '{dataset_name}' to '{output_path}'")
        else:
            print(f"Dataset '{dataset_name}' is None and was not saved.")








########### calculate rolling average for each day


def download_cmip6_netcdf_daily(start_year, end_year, model, scenario, variables, output_folder, bounds=None, num_workers=None, default_crs=None, target_dataset_crs=None):
    start_year = int(start_year)
    end_year = int(end_year)
    num_workers = int(num_workers) if num_workers is not None else 1
    
    valid_models = ["ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CESM2-WACCM",
                    "CESM2", "CMCC-CM2-SR5", "CMCC-ESM2", "CNRM-CM6-1",
                    "CNRM-ESM2-1", "EC-Earth3-Veg-LR", "EC-Earth3", "FGOALS-g3",
                    "GFDL-CM4", "GFDL-CM4_gr1", "GFDL-ESM4", "GISS-E2-1-G",
                    "HadGEM3-GC31-LL", "HadGEM3-GC31-MM", "IITM-ESM", "INM-CM4-8",
                    "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "KIOST-ESM",
                    "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
                    "MRI-ESM2-0", "NESM3", "NorESM2-LM", "NorESM2-MM", "TaiESM1",
                    "UKESM1-0-LL"]
    
    if model not in valid_models:
        raise ValueError(f"Invalid model '{model}'. Choose from {valid_models}.")

    valid_scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
    if scenario not in valid_scenarios:
        raise ValueError(f"Invalid scenario '{scenario}'. Choose from {valid_scenarios}.")

    fs = s3fs.S3FileSystem(anon=True)

    ensemble_mapping = {
        "CESM2": "r4i1p1f1",
        "CNRM-CM6-1": "r1i1p1f2",
        "GISS-E2-1-G": "r1i1p1f2",
        "MIROC-ES2L": "r1i1p1f2",
        "UKESM1-0-LL": "r1i1p1f2",
        "FGOALS-g3": "r3i1p1f1",
        "HadGEM3-GC31-LL": "r1i1p1f3",
        "HadGEM3-GC31-MM": "r1i1p1f3"
    }
    ensemble = ensemble_mapping.get(model, "r1i1p1f1")

    def get_base_path(variable):
        return f's3://nex-gddp-cmip6/NEX-GDDP-CMIP6/{model}/{scenario}/{ensemble}/{variable}/'

    def process_file(file_path, variable, fs, default_crs, target_dataset_crs, output_path):
        try:
            with fs.open(file_path, mode='rb') as f:
                ds = xr.open_dataset(f, engine='h5netcdf')
                ds.load()
                
                # Check if CRS is set and assign default_crs if not set
                if not ds.rio.crs:
                    if default_crs:
                        ds.rio.write_crs(default_crs, inplace=True)
                
                # Re-project to target_dataset_crs if specified and different
                if target_dataset_crs and ds.rio.crs != target_dataset_crs:
                    ds = ds.rio.reproject(target_dataset_crs)
                
                if bounds:
                    ds = ds.rio.clip_box(*bounds)

                # Save the dataset to the specified output path
                ds.to_netcdf(output_path)
                print(f"Saved dataset to {output_path}")
                
            return ds

        except Exception as e:
            print(f"Failed to process file {file_path}: {e}")
            return None

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        for variable in variables:
            base_path = get_base_path(variable)
            futures = []
            for year in range(start_year, end_year + 1):
                file_path = f"{base_path}{variable}_day_{model}_{scenario}_{ensemble}_gn_{year}.nc"
                output_path = os.path.join(output_folder, scenario, model, variable, f"{variable}_day_{model}_{scenario}_{ensemble}_gn_{year}.nc")
                os.makedirs(os.path.dirname(output_path), exist_ok=True)
                if "_v1.1" not in file_path:
                    future = executor.submit(process_file, file_path, variable, fs, default_crs, target_dataset_crs, output_path)
                    futures.append(future)

            for future in as_completed(futures):
                future.result()  # Wait for all tasks to complete

    print(f"All datasets for model {model} and scenario {scenario} processed successfully.")


#######
## 

def download_multiple_scenarios_daily(start_year, end_year, model, variables, scenarios, output_folder, bounds=None, num_workers=4, log_func=None, default_crs=None, target_dataset_crs=None):
    for scenario in scenarios:
        key = f"{model}_{scenario}"
        message = f"Starting processing for: {key}"
        if log_func:
            log_func(message)
        else:
            print(message)
        try:
            download_cmip6_netcdf_daily(start_year, end_year, model, scenario, variables, output_folder, bounds, num_workers, default_crs, target_dataset_crs)
            message = f"Data processed for: {key}"
            if log_func:
                log_func(message)
            else:
                print(message)
        except Exception as e:
            message = f"Exception while processing {key}: {str(e)}"
            if log_func:
                log_func(message)
            else:
                print(message)



##########
### 

def download_multiple_models_daily(start_year, end_year, models, variables, scenarios, output_folder, bounds=None, num_workers=4, log_func=None, default_crs=None, target_dataset_crs=None):
    for model in models:
        message = f"Starting processing for model: {model}"
        if log_func:
            log_func(message)
        else:
            print(message)
        
        download_multiple_scenarios_daily(start_year, end_year, model, variables, scenarios, output_folder, bounds, num_workers, log_func, default_crs, target_dataset_crs)
        
        message = f"Finished processing for model: {model}"
        if log_func:
            log_func(message)
        else:
            print(message)




#####
## aggregate daily data 


def aggregate_daily_data(input_folder, output_file, target_year):
    # List all the NetCDF files in the input folder
    file_list = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.nc')]
    
    if not file_list:
        print("No NetCDF files found in the input folder.")
        return
    
    # Open all the datasets
    datasets = [xr.open_dataset(file) for file in file_list]
    
    # Concatenate along the time dimension
    combined = xr.concat(datasets, dim="time")
    
    # Compute the mean for each day of the year
    daily_mean = combined.groupby("time.dayofyear").mean(dim="time")
    
    # Create a new time coordinate for a single year with 365 days
    single_year = xr.cftime_range(start=f"{target_year}-01-01", end=f"{target_year}-12-31", calendar="noleap")
    daily_mean = daily_mean.assign_coords(time=single_year)
    
    # Save the aggregated dataset to a NetCDF file
    daily_mean.to_netcdf(output_file)
    print(f"Aggregated data saved to {output_file}")

