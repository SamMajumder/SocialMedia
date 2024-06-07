# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 07:25:10 2024

@author: Dr. Sambadi Majumder
"""

from utils import *


def main():
    input_folder_name = input("What is the name of the folder where the data is downloaded?: ")
    scenarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
    
    model_names_input = input("What models were downloaded? Please write each model name separated with a comma: ")
    model_names = [model.strip() for model in model_names_input.split(",")]
    
    variables_input = input("What climate variables were downloaded? Please write them separated with a comma: ")
    variables = [variable.strip() for variable in variables_input.split(",")]
    
    target_year = int(input("Please specify the target year for the aggregated data: "))
    
    base_input_folder = os.path.join("..", input_folder_name)
    
    output_folder_name = input("Please state the name of the output folder you'd like: ")
    base_output_folder = os.path.join("..", output_folder_name)
    
    # Create the output folder if it does not exist
    if not os.path.exists(base_output_folder):
        os.makedirs(base_output_folder)
    
    for scenario in scenarios:
        for model in model_names:
            for variable in variables:
                input_folder = os.path.join(base_input_folder, scenario, model, variable)
                output_subfolder = os.path.join(base_output_folder, scenario, model, variable)
                output_file = os.path.join(output_subfolder, f"{variable}_aggregated_daily_{target_year}.nc")
                
                if os.path.exists(input_folder):
                    if not os.path.exists(output_subfolder):
                        os.makedirs(output_subfolder)
                    print(f"Aggregating data for {scenario}/{model}/{variable}...")
                    aggregate_daily_data(input_folder, output_file, target_year)
                else:
                    print(f"Directory {input_folder} does not exist. Skipping...")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"An unexpected error occurred: {e}")





