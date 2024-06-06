# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 08:58:13 2024

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
    
    try:
        # Extract bounds
        bounds = extract_bounds_osmnx(region_of_interest)
    
        # Ask the user for a default CRS and target CRS
        default_crs = input("Enter the default CRS to use if CRS is not present in the dataset (e.g., 'EPSG:4326'): ")
        target_dataset_crs = input("Enter the target CRS to reproject the dataset (e.g., 'EPSG:4326'): ")
        
        # Ask the user for an output directory path
        output_folder = input("Enter the name of the output folder: ")
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
        
        # Retrieve the data
        download_multiple_models_daily(
            start_year=start_year,
            end_year=end_year,
            models=climate_models,
            variables=climate_variables,
            scenarios=scenarios,
            output_folder=output_folder,
            bounds=bounds,
            num_workers=4,
            log_func=log_message,
            default_crs=default_crs,
            target_dataset_crs=target_dataset_crs
        )  
    except Exception as e:
        print(f"An error occurred during data processing: {e}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
