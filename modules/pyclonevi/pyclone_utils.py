#!/usr/bin/env python3

import sys
import pandas as pd
import argparse

# Define your functions here
# def add(a, b):
#     """Add two numbers."""
#     return a + b


#input data preprocessing
def create_pyclone_input(input_data, patient_id,output_data):
    df = pd.read_csv(input_data, sep = '\t')
    df['normal_cn'] = 2
    df['mutation_id'] = df.apply(lambda row: f"{patient_id}:{row['chr'][3:]}:{row['from']}:{row['alt']}", axis=1)
    df[['major_cn', 'minor_cn']] = df['karyotype'].str.split(':', expand=True)
    df = df.drop(columns=['karyotype'])
    df = df.rename(columns={'Indiv': 'sample_id', 'NR': 'ref_counts', 'NV': 'alt_counts', 'purity': 'tumour_content'})
    df = df.loc[:, ['mutation_id', 'sample_id', 'ref_counts', 'alt_counts', 'normal_cn', 'major_cn','minor_cn','tumour_content']]
    df.to_csv(output_data, sep="\t",index=False, header=True)
    return df

# def divide(a, b):
#     """Divide two numbers."""
#     if b == 0:
#         raise ValueError("Cannot divide by zero!")
#     return a / b

def main():
    # Create the top-level parser
    parser = argparse.ArgumentParser(description='Calculator script')
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Create sub-parser for the "create_pyclone_input" command
    parser_create_pyclone_input = subparsers.add_parser('create_pyclone_input', help='Add two numbers')
    parser_create_pyclone_input.add_argument('input_data', help='Joint table path')
    parser_create_pyclone_input.add_argument('patient_id', help='Patient identifier')
    parser_create_pyclone_input.add_argument('output_data', help='Output path for new table')



    # # Create sub-parser for the "divide" command
    # parser_divide = subparsers.add_parser('divide', help='Divide two numbers')
    # parser_divide.add_argument('a', type=float, help='First number')
    # parser_divide.add_argument('b', type=float, help='Second number')

    # Parse the arguments
    args = parser.parse_args()

    # Call the appropriate function based on the command
    if args.command == 'create_pyclone_input':
        result = add(args.input_data, args.patient_id, args.output_data)
    # elif args.command == 'subtract':
    #     result = subtract(args.a, args.b)
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
