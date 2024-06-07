# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 08:05:10 2024

@author: Dr. Sambadi Majumder
"""

from utils import * 

def log_message(message):
    print(message)

def main():
    valid_models = [
        "ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CESM2-WACCM",
        "CESM2", "CMCC-CM2-SR5", "CMCC-ESM2", "CNRM-CM6-1", "CNRM-ESM2-1",
        "EC-Earth3-Veg-LR", "EC-Earth3", "FGOALS-g3", "GFDL-CM4", "GFDL-CM4_gr1",
        "GFDL-ESM4", "GISS-E2-1-G", "HadGEM3-GC31-LL", "HadGEM3-GC31-MM",
        "IITM-ESM", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G",
        "KIOST-ESM", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
        "MRI-ESM2-0", "NESM3", "NorESM2-LM", "NorESM2-MM", "TaiESM1",
        "UKESM1-0-LL"
    ]

    valid_climate_variables = [
        "tas", "tasmax", "tasmin", "pr", "hurs", "huss", "rlds", "rsds", "sfcWind"
    ]
    
    region_of_interest = input("Please enter the region you are interested in (country or continent): ")
    
    while True:
        climate_model_input = input("Please enter climate model/s (separate multiple models with commas): ")
        climate_models = [model.strip() for model in climate_model_input.split(",") if model.strip() in valid_models]
        if climate_models:
            break
        else:
            print(f"Invalid input. Please enter valid climate models from: {valid_models}")
    
    while True:
        climate_variable_input = input("Please enter climate variable/s (separate multiple variables with commas): ")
        climate_variables = [variable.strip() for variable in climate_variable_input.split(",") if variable.strip() in valid_climate_variables]
        if climate_variables:
            break
        else:
            print(f"Invalid input. Please enter valid climate variables from: {valid_climate_variables}")
    
    while True:
        try:
            target_year = int(input("Please specify a target year between 2023 and 2090: "))
            if 2023 <= target_year <= 2090:
                break
            else:
                print("Please enter a year between 2023 and 2090.")
        except ValueError:
            print("Invalid input for target year. Please enter a valid year.")
    
    # Calculate start and end year
    start_year = target_year - 9
    end_year = target_year + 10
    scenarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
    
    default_crs = input("Enter the default CRS (e.g., 'epsg:4326') to be assigned if not present (leave empty if not needed): ")
    target_dataset_crs = input("Enter the target CRS (e.g., 'epsg:4326') to reproject the dataset (leave empty if not needed): ")

    try:
        # Extract bounds
        bounds = extract_bounds_osmnx(region_of_interest)
    
        # Retrieve the data
        results = process_multiple_models_monthly(
            start_year=start_year,
            end_year=end_year,
            models=climate_models,
            variables=climate_variables,
            scenarios=scenarios,
            bounds=bounds,
            num_workers=4,
            log_func=log_message,
            default_crs=default_crs if default_crs else None,
            target_dataset_crs=target_dataset_crs if target_dataset_crs else None
        )  
    except Exception as e:
        print(f"An error occurred during data processing: {e}")
        return
    
    # Ask the user for an output directory path
    output_folder = input("Enter the path to the output folder: ")
    output_folder = os.path.join("..", output_folder)
    try:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
            print(f"Created the directory: {output_folder}")
        else:
            print(f"Directory already exists: {output_folder}")
    except Exception as e:
        print(f"Failed to create or access directory {output_folder}: {e}")
        return
    
    try:
        export_datasets_to_netcdf(results, output_folder)
    except Exception as e:
        print(f"Failed to export datasets: {e}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

