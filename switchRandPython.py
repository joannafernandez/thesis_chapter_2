#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 22:40:12 2024

@author: patricfernandez
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd

import os
import pandas as pd

def replace_dot_in_column(directory_path, output_path):
    # Get a list of text files in the directory
    file_list = [f for f in os.listdir(directory_path) if f.endswith('.txt')]

    for file_name in file_list:
        # Construct the full path for each file
        file_path = os.path.join(directory_path, file_name)

        # Read the file into a DataFrame
        df = pd.read_csv(file_path, sep='\t', header=None)

        # Replace "." with "a" in the first column
        df[0] = df[0].str.replace('.', 'a')

        # Save the modified DataFrame to a new file
        output_file_path = os.path.join(output_path, file_name)
        df.to_csv(output_file_path, sep='\t', header=None, index=None)

# Example usage
directory_path = "/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit"
output_path = "/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone"
replace_dot_in_column(directory_path, output_path)

directory_path = "/Users/patricfernandez/Documents/python/rna_seq_rep_1_and_2_no_thia/"
output_path = "/Users/patricfernandez/Documents/python/rna_seq_rep_1_and_2_no_thia/herewegoagain/"
replace_dot_in_column(directory_path, output_path)


directory_path = "/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/"
output_path = "/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/firstattempt/"
replace_dot_in_column(directory_path, output_path)
#%%
#check 
def read_and_process_text_files(directory_path, output_path):
    # Get a list of text files in the directory
    file_list = [f for f in os.listdir(directory_path) if f.endswith('.txt')]

    # Initialize an empty dictionary to store data frames
    df_dict = {}

    # Loop through each file and read it into a data frame
    for file in file_list:
        # Extract the file name without extension to use as a data frame name
        df_name = os.path.splitext(file)[0]

        # Read the file into a data frame without header and create default column names
        file_path = os.path.join(directory_path, file)
        df = pd.read_csv(file_path, sep='\t', header=None, names=[f"col{i}" for i in range(1, 4)])
 


        # Save the modified DataFrame as a new text file
        output_path = os.path.join(output_directory, f"{df_name}_modified.txt")
        df.to_csv(output_path, sep='\t', header=None, index=False)

        # Add the data frame to the dictionary with the extracted name
        df_dict[df_name] = df

    # Return the dictionary of data frames
    return df_dict

result_dict2 = read_and_process_text_files(directory_path, output_path)

directory_path = "/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone"


#%%

def read_and_process_text_files(directory_path):
    # Get a list of text files in the directory
    file_list = [f for f in os.listdir(directory_path) if f.endswith('.txt')]
    
    # Initialize an empty dictionary to store data frames
    df_dict = {}
    
    # Loop through each file and read it into a data frame
    for file in file_list:
        # Extract the file name without extension to use as a data frame name
        df_name = os.path.splitext(file)[0]
        
        # Read the file into a data frame without header and create default column names
        file_path = os.path.join(directory_path, file)
        df = pd.read_csv(file_path, sep='\t', header=None, names=[f"col{i}" for i in range(1, 4)])
        
        # Add the data frame to the dictionary with the extracted name
        df_dict[df_name] = df
    
    # Find the common set of rows among all data frames
    common_rows = set.intersection(*[set(df.index) for df in df_dict.values()])
    
    # Filter each data frame to keep only the common rows
    for df_name, df in df_dict.items():
        df_dict[df_name] = df.loc[list(common_rows)]
    
    return df_dict

# Replace 'directory_path' and 'output_directory' with your actual paths
directory_path = '/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone'
output_directory = '/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/pleaseplease'

result_dict2 = read_and_process_text_files(directory_path)

# Optionally, save the modified data frames to new text files
for df_name, df in result_dict2.items():
    output_path = os.path.join(output_directory, f"{df_name}_modified.txt")
    df.to_csv(output_path, sep='\t', header=None, index=False)
    
#%%


import os
import pandas as pd


def read_and_process_text_files(directory_path, output_directory):
    # Get a list of text files in the directory
    file_list = [f for f in os.listdir(directory_path) if f.endswith('.txt')]

    # Initialize an empty dictionary to store data frames
    df_dict = {}

    # Loop through each file and read it into a data frame
    for file in file_list:
        # Extract the file name without extension to use as a data frame name
        df_name = os.path.splitext(file)[0]

        # Read the file into a data frame without header and create default column names
        file_path = os.path.join(directory_path, file)
        df = pd.read_csv(file_path, sep='\t', header=None, names=[f"col{i}" for i in range(1, 4)])

        # Add the data frame to the dictionary with the extracted name
        df_dict[df_name] = df

    # Find the common set of rows among all data frames
    common_rows = set.intersection(*[set(df.iloc[:, 0]) for df in df_dict.values()])

    # Filter each data frame to keep only the common rows
    for df_name, df in df_dict.items():
        df_dict[df_name] = df[df.iloc[:, 0].isin(common_rows)]

    # Ensure the first column is in the same order across all data frames
    reference_order = df_dict[list(df_dict.keys())[0]].iloc[:, 0]
    for df_name, df in df_dict.items():
        # Drop the existing col1 before resetting the index
        df_dict[df_name] = df.set_index(df.iloc[:, 0]).reindex(reference_order).drop(columns='col1').reset_index()

    # Drop the entire 'col3' column
    for df_name, df in df_dict.items():
        df_dict[df_name] = df.drop(columns='col2')

    # Save the modified data frames to new text files
    for df_name, df in df_dict.items():
        output_path = os.path.join(output_directory, f"{df_name}_modified.txt")
        df.to_csv(output_path, sep='\t', header=None, index=False)

    return df_dict

# Replace 'directory_path' and 'output_directory' with your actual paths
directory_path = '/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone'
output_directory = '/Users/patricfernandez/Documents/python/RNA_seq_rep1_edit/thisone/pleaseplease'



directory_path = "/Users/patricfernandez/Documents/python/rna_seq_rep_1_and_2_no_thia/herewegoagain/"
output_directory = "/Users/patricfernandez/Documents/python/rna_seq_rep_1_and_2_no_thia/herewegoagain/okaynowdrop/"


result_dict2 = read_and_process_text_files(directory_path, output_directory)


#this is for the combined thiamine and no thiamine

directory_path = "/Users/patricfernandez/Documents/python/rna_seq_rep_1_and_2_no_thia/herewegoagain/allofthem/"
output_directory = "/Users/patricfernandez/Documents/python/rna_seq_rep_1_and_2_no_thia/herewegoagain/allofthem/theoutput/"



result_dict2 = read_and_process_text_files(directory_path, output_directory)
df1 = result_dict2["RZ316_1_untreated"]


directory_path = "/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/firstattempt/"
output_directory = "/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/firstattempt/tryingoneofthefuncs/"
result_dict2 = read_and_process_text_files(directory_path, output_directory)
#%%

#now an issue arises with the fraction method that apparently DEseq doesn't like non-integers
#wish someone told me that fucking earlier. so chatgpt has written this for me instead:


def process_files(directory):
    # List all files in the directory
    files = [file for file in os.listdir(directory) if file.endswith('.txt')]

    for file in files:
        file_path = os.path.join(directory, file)

        # Read the file into a pandas DataFrame
        df = pd.read_csv(file_path, sep='\t', header=None, names=['ID', 'Count'])

        # Convert the 'Count' column to integers (rounding if necessary)
        df['Count'] = df['Count'].astype(int)

        # Overwrite the original file with the updated DataFrame
        df.to_csv(file_path, sep='\t', index=False, header=False)

# Replace 'your_directory_path' with the actual directory path
directory_path = '/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/firstattempt/rename_for_R/wt_and_double'
directory_path = '/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/firstattempt/rename_for_R/wt_and_dbl'
directory_path = '/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/firstattempt/rename_for_R/wt_and_sen'

directory_path = '/Users/patricfernandez/Documents/python/htseqoutput_29.1.24/firstattempt/rename_for_R/with_thia_all/allwithout_forprocessing/'
process_files(directory_path)


