#!/usr/bin/env python3

import sys
import pandas as pd
import argparse
import h5py
import numpy as np

#input data preprocessing
def create_pyclone_input(input_data, patient_id,output_data):
    df = pd.read_csv(input_data, sep = '\t',header=0)
    df['normal_cn'] = 2
    df['patient_id'] = patient_id
    df['mutation_id'] = df.apply(lambda row: f"{patient_id}:{row['chr'][3:]}:{row['from']}:{row['alt']}", axis=1)
    df[['major_cn', 'minor_cn']] = df['karyotype'].str.split(':', expand=True)
    df = df.drop(columns=['karyotype'])
    df = df.rename(columns={'Indiv': 'sample_id', 'NR': 'ref_counts', 'NV': 'alt_counts', 'purity': 'tumour_content'})


    column_names = ['mutation_id', 'patient_id','sample_id', 'ref_counts', 'alt_counts', 'normal_cn', 'major_cn','minor_cn','tumour_content','driver_label','is_driver']
    column_names_red = list(set(column_names) & set(df.columns))

    # df = df.loc[:, ['mutation_id', 'patient_id','sample_id', 'ref_counts', 'alt_counts', 'normal_cn', 'major_cn','minor_cn','tumour_content','driver_label','is_driver']]
    # column_names = ['mutation_id', 'patient_id','sample_id', 'ref_counts', 'alt_counts', 'normal_cn', 'major_cn','minor_cn','tumour_content','driver_label','is_driver']

    df = df.loc[:, column_names_red]
    df.to_csv(output_data, sep="\t",index=False, header=df.columns)
    print(df)

def update_hdf5_with_csv(existing_file_path: str, pyclone_tsv: str) -> None:
    """
    Reads data from a CSV file, processes it, and appends it to an existing HDF5 file.

    Parameters:
    existing_file_path (str): Path to the existing HDF5 file.
    csv_file_path (str): Path to the CSV file to be read and processed.
    """
    # Read the CSV file into a pandas DataFrame
    original_data = pd.read_csv(pyclone_tsv, sep='\t')
    original_data["VAF"] = original_data["alt_counts"] / (original_data["alt_counts"] + original_data["ref_counts"])
    
    VAF = original_data["VAF"].to_numpy()
    mutation_id = original_data["mutation_id"].to_numpy()
    sample_id = original_data["sample_id"].to_numpy()
    
    # Open the existing HDF5 file in append mode
    with h5py.File(existing_file_path, 'a') as existing_file:
        # Create a group for the original data if it doesn't already exist
        if 'original_data' not in existing_file:
            group = existing_file.create_group('original_data')
        else:
            group = existing_file['original_data']
        
        # Create datasets within the group for VAF, mutation_id, and sample_id
        group.create_dataset("VAF", data=VAF)
        group.create_dataset("mutation_id", data=mutation_id)
        group.create_dataset("sample_id", data=sample_id)


def main():
    # Create the top-level parser
    parser = argparse.ArgumentParser(description='Calculator script')
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Create sub-parser for the "create_pyclone_input" command
    parser_create_pyclone_input = subparsers.add_parser('create_pyclone_input', help='Add two numbers')
    parser_create_pyclone_input.add_argument('input_data', help='Joint table path')
    parser_create_pyclone_input.add_argument('patient_id', help='Patient identifier')
    parser_create_pyclone_input.add_argument('output_data', help='Output path for new table')

    # Create sub-parser for the "update_hdf5_with_csv" command
    parser_update_hdf5_with_csv = subparsers.add_parser('update_hd5_file_with_csv')
    parser_update_hdf5_with_csv.add_argument('existing_file_path', help='Path to H5 file')
    parser_update_hdf5_with_csv.add_argument('pyclone_tsv', help='Pyclone input file')

    # Parse the arguments
    args = parser.parse_args()

    # Call the appropriate function based on the command
    if args.command == 'create_pyclone_input':
        result = create_pyclone_input(args.input_data, args.patient_id, args.output_data)
    elif args.command == 'update_hd5_file_with_csv':
        result = update_hdf5_with_csv(args.existing_file_path, args.pyclone_tsv)
    # elif args.command == 'multiply':
    #     result = multiply(args.a, args.b)
    # elif args.command == 'divide':
    #     result = divide(args.a, args.b)
    else:
        parser.print_help()
        return

    # Print the result
    print(f"The result is: {result}")

if __name__ == '__main__':
    main()
